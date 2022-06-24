'''
Collection of routines related to Karhunen-Loeve modes.

The theory is based on the paper
"Optimal bases for wave-front simulation and
reconstruction on annular apertures", Robert C. Cannon, 1996, JOSAA, 13, 4

The present implementation is based on the IDL package of R. Cannon
(wavefront modelling and reconstruction).
A closely similar implementation can also be find in Yorick in the YAO package.

Usage
-----
Main routine is 'make_kl' to generate KL basis of dimension ``[dim, dim, nmax]``.

For Kolmogorov statistics, e.g. ::

    kl, _, _, _ = make_kl(150, 128, ri = 0.2, stf='kolmogorov')

.. warning::

    make_kl with von Karman stf fails. It has been implemented but KL generation failed in 'while loop' of gkl_fcom...

.. codeauthor:: Gilles Orban de Xivry (ULieve)
:date: November 2017
'''

import numpy as np
import scipy
from scipy.ndimage.interpolation import map_coordinates


def rebin(a, newshape):
    '''Rebin an array to a new shape.
    See scipy cookbook.
    It is intended to be similar to 'rebin' of IDL
    '''
    assert len(a.shape) == len(newshape)

    slices = [slice(0, old, float(old) / new)
              for old, new in zip(a.shape, newshape)]
    coordinates = np.mgrid[slices]
    # choose the biggest smaller integer index
    indices = coordinates.astype('i')
    return a[tuple(indices)]


def stf_kolmogorov(r):
    '''
    Kolmogorov structure function with r = (D/r0)
    '''
    return 6.8839 * r**(5. / 3)


def stf_vonKarman_yao(r, L):
    return 6.88 * r**(5. / 3) * (1 - 1.485 * (r / L)**(1 / 3.) +
                                 5.383 * (r / L)**2 -
                                 6.281 * (r / L)**(7 / 3.))


def stf_vonKarman(r, L0):
    '''
        von Karman structure function with r = (D / r0)
        L0 is in unit of telescope diameter, typically a few (3; or 20m)
    '''
    r0 = 1
    D_vk = (0.17253 * (L0 / (r0)) ** (5. / 3.)
            * (1 - 2 * np.pi ** (5. / 6.) * ((r) / L0) ** (5. / 6.)
               / scipy.special.gamma(5. / 6.)
               * scipy.special.kv(5. / 6., (2 * np.pi * r) / L0)))
    return D_vk


def gkl_radii(ri, nr):
    '''
    Generate n points evenly spaced in r^2 between r_i^2 and 1

    Parameters
    -----------
        nr : int
            number of resolution elements
        ri : float
            radius of central obscuration radius (normalized; <1)

    Returns
    -------
        r: 1d-array, float
            grid of point to calculate the appropriate kernels.
            correspond to 'sqrt(s)' wrt the paper.
    '''
    d = (1 - ri**2) / nr
    # r2 = ri**2 + d / 2. + d * np.arange(nr)
    r2 = ri**2 + d * np.arange(nr) + d / 16

    return np.sqrt(r2)


def gkl_kernel(ri, nr, rad, stfunc='kolmogorov', outerscale=None):
    '''
    Calculation of the kernel L^p

    The kernel constructed here should be simply a discretization
    of the continuous kernel.
    It needs rescaling before it is treated as a matrix for finding
    the eigen-values

    Parameters
    ----------
        ri : float
            radius of central obscuration radius (normalized by D/2; <1)
        nr : int
            number of resolution elements
        rad : 1d-array
            grid of points where the kernel is evaluated
        stfunc : string
            string tag of the structure function on which the kernel are
            computed
        outerscale : float
            in unit of telescope diameter.
            Outer-scale for von Karman structure function.

    Returns
    -------
        L^p : ndarray
            kernel L^p of dimension (nr, nr, nth), where nth is the
            azimuthal discretization (5 * nr)
    '''
    oversampling = 5
    # ovsersampling for fourier calculation
    nth = oversampling * nr
    kernel = np.zeros((nr, nr, nth))

    # 1/2 because Lp and not Kp
    fnorm = 1. / 2. * (-1) / (2 * np.pi * (1 - ri**2))

    for i in range(nr):
        for j in range(i + 1):
            radius = 0.5 * np.sqrt(rad[i]**2 + rad[j]**2 -
                                   2 * rad[i] * rad[j] *
                                   np.cos(np.arange(nth) * 2 * np.pi / nth))
            if (stfunc == 'kolmogorov') or (stfunc == 'kolstf'):
                sf = stf_kolmogorov(radius)
            elif (stfunc == 'vonKarman') or (stfunc == 'karman') or \
                    (stfunc == 'vk'):
                assert outerscale is not None
                sf = stf_vonKarman(radius, outerscale)
            else:
                raise AttributeError("Structure function not implemented")

            # tmp = fnorm * (2 * np.pi / nth) * (fft.ifft(sf) + fft.fft(sf)) / 2
            # fft.dct(sf, type=3) # fft.ifft(sf) # ift(sf, 1)
            tmp = fnorm * (2 * np.pi / nth) * np.fft.fft(sf, axis=0)

            # Kernel is symmetric
            kernel[i, j, :] = tmp
            kernel[j, i, :] = tmp

    return kernel


def piston_orth(nr):
    '''
    Unitary matrix used to filter out piston term.
    Eq. 19 in Cannon 1996.

    Parameters
    ----------
    nr : int
        number of resolution elements

    Returns
    -------
    U : 2d array
        Unitary matrix
    '''

    s = np.zeros((nr, nr))
    for j in range(nr - 1):
        rnm = 1. / np.sqrt((j + 1) * (j + 2))
        s[0:j + 1, j] = rnm
        s[j + 1, j] = (-1) * (j + 1) * rnm
    rnm = 1. / np.sqrt(nr)
    s[:, nr - 1] = rnm
    return s


def gkl_fcom(ri, kernels, nfunc, verbose=False):
    '''
    Computation of the radial eigenvectors of the KL basis.

    Obtained by taking the eigenvectors from the matrix L^p.
    The final function corresponds to the 'nfunc' largest eigenvalues.
    See eq. 16-18 in Cannon 1996.

    Parameters
    ----------
        ri : float
            radius of central obscuration radius (normalized by D/2; <1)
        kernels : ndarray
            kernel L^p of dimension (nr, nr, nth), where nth is the
            azimuthal discretization (5 * nr)
        nfunc : int
            number of final KL functions
    Returns
    -------
        evals : 1darray
            eigenvalues
        nord : int
            resulting number of azimuthal orders
        npo : int
        ord : 1darray
        rabas : ndarray
            radial eigenvectors of the KL basis
    '''
    kers = np.copy(kernels)
    s = np.shape(kers)
    nr = s[0]
    nt = s[2]
    nxt = 0
    fktom = (1. - ri**2) / nr
    fevtos = np.sqrt(2 * nr)
    evs = np.empty((nr, nt))

    # * Zero-order is a special case *
    # see eq 19 and 20 of Cannon's paper
    zom = kers[:, :, 0]
    s = piston_orth(nr).T
    # b1 = np.dot(np.dot(s.T, zom), s)[0:nr - 1, 0:nr - 1]
    # b1 = (np.inner(s, np.inner(zom, s).T))[0:nr - 1, 0:nr - 1]
    b1 = np.dot(np.dot(s, zom), s.T)[0:nr - 1, 0:nr - 1]

    # since matrix is symmetric, can use eigh instead of svd.
    # should not make any difference however...
    newev, v0 = np.linalg.eigh(fktom * b1)
    v1 = np.zeros((nr, nr))
    v1[0:nr - 1, 0:nr - 1] = v0.T
    # v1[nr - 1, nr - 1] = 1.0
    v1[nr - 1, nr - 1] = 1.0
    vs = np.dot(v1, s)
    newev = np.append(newev, 0)
    evs[:, nxt] = newev
    kers[:, :, nxt] = np.sqrt(nr) * vs.T

    # * Other orders - more straightforward
    nxt = 1
    while True:
        newev, vs = np.linalg.eigh(fktom * kers[:, :, nxt])
        evs[:, nxt] = newev
        # kers[:, :, nxt] = np.sqrt(2 * nr) * vs
        kers[:, :, nxt] = fevtos * vs
        if verbose:
            print('{0:.4f}'.format(nxt))
        mxn = np.max(newev)
        egtmxn = np.array(np.floor(evs[:, 0:nxt + 1] > mxn), dtype='int')
        nxt = nxt + 1
        if ((2 * np.sum(egtmxn) - np.sum(egtmxn[:, 0])) >= nfunc):
            break
    nus = nxt - 1

    # * The rest is about sorting and selecting the N functions with the
    # highest eigenvalues *
    kers = kers[:, :, 0:nus]
    evs = np.reshape(evs[:, 0:nus].T, nr * nus)
    a = (np.argsort(-1 * evs))[0:nfunc]

    # every eigenvalue occurs twice except those for the zeroth order mode.
    # this could be done without the loops,
    # but it isn't the sticking point anyway...
    no = 0
    ni = 0
    # oind = np.empty(nfunc + 1)
    oind = np.zeros(nfunc + 1, dtype='int')
    while True:
        if (a[ni] < nr):
            oind[no] = a[ni]
            no = no + 1
        else:
            oind[no] = a[ni]
            oind[no + 1] = a[ni]
            no = no + 2
        ni = ni + 1
        if (no >= nfunc):
            break

    oind = oind[0:nfunc]
    tord = oind // nr
    odd = ((np.arange(nfunc) % 2) == 1)
    pio = oind % nr

    evals = evs[oind]
    oord = 2 * tord - np.array(np.floor((tord >= 1) & odd), dtype='int')
    nord = np.max(oord) + 1
    rabas = np.zeros((nr, nfunc))
    npo = np.zeros(np.max(oord) + 1)

    for i in range(nfunc):
        npo[oord[i]] = npo[oord[i]] + 1
        rabas[:, i] = kers[:, pio[i], tord[i]]

    return evals, nord, npo, oord, rabas


def gkl_azimuthal(nord, npp):
    '''
    Compute the azimuthal function of the KL basis.

    Parameters
    ----------
        nord : int
            number of azimuthal orders
        npp : int
            grid of point sampling the azimuthal coordinate
    Returns
    -------
        gklazi: ndarray
            azimuthal function of the KL basis.
    '''
    gklazi = np.zeros((1 + nord, npp))
    theta = np.arange(npp) * (2 * np.pi / npp)
    gklazi[0, :] = 1.0
    for i in range(1, nord, 2):
        # even
        gklazi[i, :] = np.cos((i // 2 + 1) * theta)
    for i in range(2, nord, 2):
        # odd
        gklazi[i, :] = np.sin((i // 2) * theta)
    return gklazi


def gkl_basis(ri=0.25, nr=40, npp=None, nfunc=500,
              stf='kolstf', outerscale=None):
    '''
    Wrapper to create the radial and azimuthal K-L functions.

    Parameters
    ----------
        ri : float
            normalized internal radius
        nr : int
            number of radial resolution elements
        np : int
            number of azimuthal resolution elements
        nfunc : int
            number of generated K-L function
        stf : string
            structure function tag describing the atmospheric statistics
    Returns
    -------
        gklbasis: dic
            dictionary containing the radial and azimuthal basis + other
            relevant information
    '''
    if npp is None:
        npp = 5 * nr

    if (nr * npp) / nfunc < 8:
        print("warning: you may need a finer radial sampling ")
        print("(ie, increased 'nr') to generate {0:1.0f} "
              "functions".format(nfunc))
    elif (nr * npp) / nfunc > 40:
        print("note, for this size basis, radial discretization on ", nr)
        print("points is finer than necessary - it should work, but you ")
        print("could take a smaller 'nr' without loss of accuracy")

    rad_basis = gkl_radii(ri, nr)

    kers = gkl_kernel(ri, nr, rad_basis, stf, outerscale)

    evals, nord, npo, oord, rabas = gkl_fcom(ri, kers, nfunc)

    azi_basis = gkl_azimuthal(nord, npp)

    gklbasis = {'nr': nr, 'np': npp, 'nfunc': nfunc, 'ri': ri, 'stfn': ' ',
                'radp': rad_basis, 'evals': evals,
                'nord': nord, 'npo': npo, 'ord': oord,
                'rabas': rabas, 'azbas': azi_basis}
    return gklbasis


def gkl_sfi(kl_basis, i):
    '''
    return the i'th function from the generalized KL basis 'bas'.
    'bas' must be generated by 'gkl_basis'
    '''
    if i > kl_basis['nfunc']:
        raise("The basis only contains {0:1.0f} "
              "functions".format(kl_basis['nfunc']))
    nr = kl_basis['nr']
    npp = kl_basis['np']
    oord = kl_basis['ord'][i]

    rad_bas = rebin(np.reshape(kl_basis['rabas'][:, i], (nr, 1)), (nr, npp))
    az_bas = rebin(np.reshape(kl_basis['azbas'][oord, :], (1, npp)), (nr, npp))

    sf = rad_bas * az_bas

    return sf


def radii(nr, npp, ri):
    '''
    Use to generate a polar coordinate system.

    Generate an nr x npp array with npp copies of the radial coordinate
    array.
    Radial coordinate span the range from r=ri to r=1 with
    successive annuli having equal areas.

    ie, the area between ri and 1 is divided into nr equal rings, and
    the points are positioned at the half-area mark on each ring;
    there are no points on the border

    see also
    --------
    polang
    '''
    # r2 = ri**2 + (np.arange(nr) + 0.5) / nr * (1 - ri**2)
    r2 = ri**2 + (np.arange(nr)) / nr * (1 - ri**2)
    rs = np.sqrt(r2)
    ra = rebin(np.reshape(rs, (nr, 1)), (nr, npp))
    # ra = rebin(rs, nr, npp)
    return ra


def polang(r):
    '''
        Generate an array with the same dimensions as r, but containing the
        azimuthal values for a polar coordinate system.
    '''
    s = np.shape(r)
    nr = s[0]
    npp = s[1]
    # phi = np.zeros((nr, npp))
    phi1 = np.arange(npp) / npp * 2.0 * np.pi

    phi = np.transpose(rebin(np.reshape(phi1, (npp, 1)), (npp, nr)))
    return phi


def set_pctr(bas, ncp=None, ncmar=None):
    '''
    call pcgeom to build the dic geom_struct with the
    right initializtions.

    Parameters
    ----------
        bas : dic
            gkl_basis dic built with the gkl_bas routine
    '''
    if ncmar is None:
        ncmar = 2
    if ncp is None:
        ncp = 128

    return pcgeom(bas['nr'], bas['np'], ncp, bas['ri'], ncmar)


def setpincs(ax, ay, px, py, ri):
    '''
    determine a set of squares for interpolating from cartesian
    to polar coordinates,
    using only those points with ri <= r <= 1
    '''
    s = np.shape(ax)
    nc = s[0]
    s = np.shape(px)
    nr = s[0]
    npp = s[1]
    # dcar = (ax[0, nc - 1] - ax[0, 0]) / (nc - 1)
    dcar = (ax[0, nc - 1] - ax[0, 0]) / (nc - 1)
    ofcar = ax[0, 0]

    rlx = (px - ofcar) / dcar
    rly = (py - ofcar) / dcar
    lx = np.array(rlx, dtype='int')
    ly = np.array(rly, dtype='int')
    # shx = rlx - lx
    # shy = rly - ly
    shx = rlx - lx - 1
    shy = rly - ly - 1
    #
    # pincx = np.array([lx.T, (lx + 1).T, (lx + 1).T, lx.T])
    # pincy = np.array([ly.T, ly.T, (ly + 1).T, (ly + 1).T])
    pincx = np.array([(lx - 1).T, (lx).T, (lx).T, (lx - 1).T])
    pincy = np.array([(ly - 1).T, (ly - 1).T, (ly).T, (ly).T])
    pincw = np.array([(1 - shx).T, (shx).T, (shx).T, (1 - shx).T]) * \
        np.array([(1 - shy).T, (1 - shy).T, shy.T, shy.T])
    axy = ax**2 + ay**2
    axyinap = ((axy >= ri**2.) & (axy <= 1.0))
    # axyinap = ((axy >= ri**2.) & (axy <= 1.0))
    # print(pincx.shape, pincy.shape, axyinap.shape, pincw.shape)
    pincw = np.squeeze(pincw * axyinap[pincx, pincy])
    pincw = pincw * rebin(np.reshape(1.0 / np.sum(pincw, 0), (1, nr, npp)),
                          (4, npp, nr))

    return pincx, pincy, pincw


def pcgeom(nr, npp, ncp, ri, ncmar):
    '''
    This routine builds a geom dic.

    Parameters
    ----------
        px, py : int
            the x, y coordinates of points in the polar arrays.
        cr, cp :
            the r, phi coordinates of points in the cartesian grids.
        ncmar :
            allows the possibility that there is a margin of ncmar points in the
            cartesian arrays outside the region of interest.
    '''
    nused = ncp - 2 * ncmar
    ff = 0.5 * nused
    hw = float(ncp - 1) / 2
    # hw = float(ncp) / 2

    r = radii(nr, npp, ri)
    p = polang(r)

    px0 = r * np.cos(p)
    py0 = r * np.sin(p)
    px = ff * px0 + hw
    py = ff * py0 + hw
    ax = np.reshape(np.arange(ncp * ncp), (ncp, ncp)) % ncp - 0.5 * (ncp - 1)
    # ax = ax / (0.5 * nused - 0.5)
    ax = ax / (0.5 * nused)
    ay = np.transpose(ax)
    pincx, pincy, pincw = setpincs(ax, ay, px0, py0, ri)
    dpi = 2 * np.pi
    cr2 = (ax**2 + ay**2)
    # ap = (cr2 > ri**2) & (cr2 < 1.)
    ap = (cr2 >= ri**2) & (cr2 <= 1.)
    # ap = np.clip(cr2, ri**2 + 1e-3, 0.999)
    # cr = (cr2 - ri**2) / (1 - ri**2) * nr - 0.5
    cr = (cr2 - ri**2) / (1 - ri**2) * nr  # - 0.5
    cp = (np.arctan2(ay, ax) + dpi) % dpi
    cp = (npp / dpi) * cp

    # cr = np.clip(cr, 1.e-3, nr - 1.001)
    # cp = np.clip(cp, 1.e-3, npp - 1.001)
    cr = np.clip(cr, 1e-3, nr - 1.001)  # - 1.00)
    cp = np.clip(cp, 1e-3, npp - 1.001)  # - 1.00)

    geom = {'px': px, 'py': py, 'cr': cr, 'cp': cp,
            'pincx': pincx, 'pincy': pincy, 'pincw': pincw,
            'ap': ap, 'ncp': ncp, 'ncmar': ncmar}

    return geom


def pol2car(cpgeom, pol, mask=False):
    '''
    Polar to cartesian conversion.

    (points not in the aperture are treated as though they were at
    the first or last radial polar value...)
    '''
    # f = interp2d(cpgeom['cr'], cpgeom['cp'], pol)

    cd = map_coordinates(pol, [cpgeom['cr'], cpgeom['cp']],
                         order=1, mode='nearest')
    if mask is not False:
        cd = cd * cpgeom['ap']
    return cd


def make_kl(nmax, dim, ri=0.0, nr=40,
            stf='kolmogorov', outerscale=None, mask=True):
    '''
    Main routine to generatre a KL basis of dimension [nmax, dim, dim].

    For Kolmogorov statistics, e.g. ::

        kl, _, _, _ = make_kl(150, 128, ri = 0.2, stf='kolmogorov')

    | As a rule of thumb
    |     nr x npp = 50 x 250 is fine up to 500 functions
    |     60 x 300 for a thousand
    |     80 x 400 for three thousands.

    Parameters
    ----------
    nmax : int
        number of KL function to generate
    dim : int
        size of the KL arrays
    ri : float
        radial central obscuration normalized by D/2
    nr : int
        number of point on radius. npp (number of azimuthal pts) is  =2 pi * nr
    stf : string
        structure function tag. Default is 'kolmogorov'
    outerscale : float
        outer scale in units of telescope diameter. Releveant if von vonKarman
        stf. (not implemented yet)
    mask : bool
        pupil masking. Default is True

    Returns
    -------
    kl : ndarray (nmax, dim, dim)
        KL basis in cartesian coordinates
    varKL : 1darray
        associated variance
    pupil : 2darray
        pupil
    polar_base : dictionary
        polar base dictionary used for the KL basis computation.
        as returned by 'gkl_basis'

    SEE ALSO
    --------
    gkl_basis, set_pctr, pol2car
    '''
    if (ri == 0):
        raise ValueError("KL modes generation only works for annular apertures, ri must be > 0!")

    npp = int(2 * np.pi * nr)
    if (nr * npp) < (15 * nmax):
        print('Polar grid sampling may be insufficienttoto_basis2.fits, please'
              ' consider increasing nr')

    assert (stf == 'kolmogorov') or (stf == 'kolstf') or (stf == 'vonKarman') \
        or (stf == 'karman') or (stf == 'vk')

    if stf == 'vonKarman':
        assert outerscale is not None

    # from aotools import circle
    # if ri > 0.0:
    #     pup = circle(dim / 2, dim) - circle(dim / 2 * ri, dim)
    # else:
    #     pup = circle(dim / 2, dim)

    polar_base = gkl_basis(ri, nr, npp, nfunc=nmax,
                           stf=stf, outerscale=outerscale)

    pc1 = set_pctr(polar_base, ncp=dim, ncmar=0)

    kl = np.zeros((nmax, dim, dim))

    for i in range(nmax):
        kl[i, :, :] = pol2car(pc1, gkl_sfi(polar_base, i), mask=mask)
        # if mask is True:
        #     kl[i, :, :] *= pup

    pupil = np.array(pc1['ap'], dtype='float')
    # pupil = pup
    varKL = polar_base['evals']
    return kl, varKL, pupil, polar_base
