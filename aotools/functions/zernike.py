import numpy
from . import circle

def phaseFromZernikes(zCoeffs, size, norm="noll", rot=0):
    """
    Creates an array of the sum of zernike polynomials with specified coefficeints

    Parameters:
        zCoeffs (list): zernike Coefficients
        size (int): Diameter of returned array
        norm (string, optional): The normalisation of Zernike modes. Can be ``"noll"``, ``"p2v"`` (peak to valley), or ``"rms"``. default is ``"noll"``.
        rot (float) : Rotates the Zernike mode by rot radians around its centre. Defaults to zero.
    Returns:
        ndarray: a `size` x `size` array of summed Zernike polynomials
    """
    Zs = zernikeArray(len(zCoeffs), size, norm=norm, rot=rot)
    phase = numpy.zeros((size, size))
    for z in range(len(zCoeffs)):
        phase += Zs[z] * zCoeffs[z]

    return phase


def zernike_noll(j, N, rot=0):
    """
     Creates the Zernike polynomial with mode index j,
     where j = 1 corresponds to piston.

     Args:
        j (int): The noll j number of the zernike mode
        N (int): The diameter of the zernike more in pixels
        rot (float): Rotates the Zernike mode by rot radians around its centre. Defaults to zero.
     Returns:
        ndarray: The Zernike mode
     """

    n, m = zernIndex(j)
    return zernike_nm(n, m, N, rot)


def zernike_nm(n, m, N, rot=0):
    """
     Creates the Zernike polynomial with radial index, n, and azimuthal index, m.

     Args:
        n (int): The radial order of the zernike mode
        m (int): The azimuthal order of the zernike mode
        N (int): The diameter of the zernike more in pixels
        rot (float): Rotates the zernike by rot radians around its centre. Defaults to 0.
     Returns:
        ndarray: The Zernike mode
     """
    coords = (numpy.arange(N) - N / 2. + 0.5) / (N / 2.)
    X, Y = numpy.meshgrid(coords, coords)
    R = numpy.sqrt(X**2 + Y**2)
    theta = numpy.arctan2(Y, X)

    if m==0:
        Z = numpy.sqrt(n+1)*zernikeRadialFunc(n, 0, R)
    else:
        if m > 0: # j is even
            Z = numpy.sqrt(2*(n+1)) * zernikeRadialFunc(n, m, R) * numpy.cos((m*theta)+rot)
        else:   #i is odd
            m = abs(m)
            Z = numpy.sqrt(2*(n+1)) * zernikeRadialFunc(n, m, R) * numpy.sin((m*theta)+rot)

    # clip
    Z = Z*numpy.less_equal(R, 1.0)

    return Z*circle(N/2., N)


def zernikeRadialFunc(n, m, r):
    """
    Fucntion to calculate the Zernike radial function

    Parameters:
        n (int): Zernike radial order
        m (int): Zernike azimuthal order
        r (ndarray): 2-d array of radii from the centre the array

    Returns:
        ndarray: The Zernike radial function
    """

    R = numpy.zeros(r.shape)
    # Can cast the below to "int", n,m are always *both* either even or odd
    for i in range(0, int((n - m) / 2) + 1):

        R += numpy.array(r**(n - 2 * i) * (((-1)**(i)) *
                         numpy.math.factorial(n - i)) /
                         (numpy.math.factorial(i) *
                          numpy.math.factorial(int(0.5 * (n + m) - i)) *
                          numpy.math.factorial(int(0.5 * (n - m) - i))),
                         dtype='float')
    return R


def zernIndex(j):
    """
    Find the [n,m] list giving the radial order n and azimuthal order
    of the Zernike polynomial of Noll index j.

    Parameters:
        j (int): The Noll index for Zernike polynomials

    Returns:
        list: n, m values
    """
    n = int((-1.+numpy.sqrt(8*(j-1)+1))/2.)
    p = (j-(n*(n+1))/2.)
    k = n%2
    m = int((p+k)/2.)*2 - k

    if m!=0:
        if j%2==0:
            s=1
        else:
            s=-1
        m *= s

    return [n, m]


def zernikeArray(J, N, norm="noll", rot=0):
    """
    Creates an array of Zernike Polynomials

    Parameters:
        maxJ (int or list): Max Zernike polynomial to create, or list of zernikes J indices to create
        N (int): size of created arrays
        norm (string, optional): The normalisation of Zernike modes. Can be ``"noll"``, ``"p2v"`` (peak to valley), or ``"rms"``. default is ``"noll"``.
        rot (float): Rotates the zernike modes by rot radians around their centre. Defaults to 0.

    Returns:
        ndarray: array of Zernike Polynomials
    """
    # If list, make those Zernikes
    try:
        nJ = len(J)
        Zs = numpy.empty((nJ, N, N))
        for i in range(nJ):
            Zs[i] = zernike_noll(J[i], N, rot)

    # Else, cast to int and create up to that number
    except TypeError:

        maxJ = int(numpy.round(J))
        N = int(numpy.round(N))

        Zs = numpy.empty((maxJ, N, N))

        for j in range(1, maxJ+1):
            Zs[j-1] = zernike_noll(j, N, rot)


    if norm=="p2v":
        for z in range(len(Zs)):
            Zs[z] /= (Zs[z].max()-Zs[z].min())

    elif norm=="rms":
        for z in range(len(Zs)):
            # Norm by RMS. Remember only to include circle elements in mean
            Zs[z] /= numpy.sqrt(
                    numpy.sum(Zs[z]**2)/numpy.sum(circle(N/2., N)))



    return Zs


def makegammas(nzrad):
    """
    Make "Gamma" matrices which can be used to determine first derivative
    of Zernike matrices (Noll 1976).

    Parameters:
        nzrad: Number of Zernike radial orders to calculate Gamma matrices for

    Return:
        ndarray: Array with x, then y gamma matrices
    """
    n=[0]
    m=[0]
    tt=[1]
    trig=0

    for p in range(1,nzrad+1):
        for q in range(p+1):
            if(numpy.fmod(p-q,2)==0):
                if(q>0):
                    n.append(p)
                    m.append(q)
                    trig=not(trig)
                    tt.append(trig)
                    n.append(p)
                    m.append(q)
                    trig=not(trig)
                    tt.append(trig)
                else:
                    n.append(p)
                    m.append(q)
                    tt.append(1)
                    trig=not(trig)
    nzmax=len(n)

    #for j in range(nzmax):
        #print j+1, n[j], m[j], tt[j]

    gamx = numpy.zeros((nzmax,nzmax),"float32")
    gamy = numpy.zeros((nzmax,nzmax),"float32")

    # Gamma x
    for i in range(nzmax):
        for j in range(i+1):

            # Rule a:
            if (m[i]==0 or m[j]==0):
                gamx[i,j] = numpy.sqrt(2.0)*numpy.sqrt(float(n[i]+1)*float(n[j]+1))
            else:
                gamx[i,j] = numpy.sqrt(float(n[i]+1)*float(n[j]+1))

            # Rule b:
            if m[i]==0:
                if ((j+1) % 2) == 1:
                    gamx[i,j] = 0.0
            elif m[j]==0:
                if ((i+1) % 2) == 1:
                    gamx[i,j] = 0.0
            else:
                if ( ((i+1) % 2) != ((j+1) % 2) ):
                    gamx[i,j] = 0.0

            # Rule c:
            if abs(m[j]-m[i]) != 1:
                gamx[i,j] = 0.0

            # Rule d - all elements positive therefore already true

    # Gamma y
    for i in range(nzmax):
        for j in range(i+1):

            # Rule a:
            if (m[i]==0 or m[j]==0):
                gamy[i,j] = numpy.sqrt(2.0)*numpy.sqrt(float(n[i]+1)*float(n[j]+1))
            else:
                gamy[i,j] = numpy.sqrt(float(n[i]+1)*float(n[j]+1))

            # Rule b:
            if m[i]==0:
                if ((j+1) % 2) == 0:
                    gamy[i,j] = 0.0
            elif m[j]==0:
                if ((i+1) % 2) == 0:
                    gamy[i,j] = 0.0
            else:
                if ( ((i+1) % 2) == ((j+1) % 2) ):
                    gamy[i,j] = 0.0

            # Rule c:
            if abs(m[j]-m[i]) != 1:
                gamy[i,j] = 0.0

            # Rule d:
            if m[i]==0:
                pass    # line 1
            elif m[j]==0:
                pass    # line 1
            elif m[j]==(m[i]+1):
                if ((i+1) % 2) == 1:
                    gamy[i,j] *= -1.    # line 2
            elif m[j]==(m[i]-1):
                if ((i+1) % 2) == 0:
                    gamy[i,j] *= -1.    # line 3
            else:
                pass    # line 4

    return numpy.array([gamx,gamy])
