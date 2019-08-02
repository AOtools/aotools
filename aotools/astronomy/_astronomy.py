import numpy
import astropy.constants as c
from ftplib import FTP
from io import BytesIO
from astroquery.vizier import Vizier
from .intrinsic_colors import colors_dict, colors_columns

# Source : https://www.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
# Dictionary of flux values at the top of the atmosphere
#                 band, lamda, dLamda, m=0 flux (Jy)
FLUX_DICTIONARY = {'U': [0.36, 0.15, 1810],
                'B': [0.44, 0.22, 4260],
                'V': [0.55, 0.16, 3640],
                'R': [0.64, 0.23, 3080],
                'I': [1.0, 0.19, 2550],
                'J': [1.26, 0.16, 1600],
                'H': [1.60, 0.23, 1080],
                'K': [2.22, 0.23, 670],
                'g': [0.52, 0.14, 3730],
                'r': [0.67, 0.14, 4490],
                'i': [0.79, 0.16, 4760],
                'z': [0.91, 0.13, 4810]}


def photons_per_mag(mag, mask, pixel_scale, wvlBand, exposure_time):
    """
    Calculates the photon flux for a given aperture, star magnitude and wavelength band

    Parameters:
        mag (float): Star apparent magnitude
        mask (ndarray): 2-d pupil mask array, 1 is transparent, 0 opaque
        pixel_scale (float): size in metres of each pixel in mask
        wvlBand (float): length of wavelength band in nanometres
        exposure_time (float): Exposure time in seconds

    Returns:
        float: number of photons
    """
    # Area defined in cm, so turn m to cm
    area = mask.sum() * pixel_scale ** 2 * 100 ** 2

    photonPerSecPerAreaPerWvl = 1000 * (10**(-float(mag)/2.5))

    # Wavelength defined in Angstroms
    photonPerSecPerArea = photonPerSecPerAreaPerWvl * wvlBand*10

    photonPerSec = photonPerSecPerArea * area

    photons = float(photonPerSec * exposure_time)

    return photons


def flux_AB_mag(mag, wvl, wvlBand):
    """
    Calculates the photon flux for a given magnitude in a given wave band using AB definition

    Parameters:
        mag (float): magnitude in AB system
        wvl (float): central wavelength in any length unit
        wvlBand (float): length of wavelength band in the same unit as wvl

    Returns:
        float: flux of photons (ph/m2/s)
    """
    return 10 ** (-mag / 2.5) * 10 ** (-48.60/2.5 - 7 + 4) * 1/c.h.value * wvlBand / wvl


def flux_STMAG_mag(mag, wvl, wvlBand):
    """
    Calculates the photon flux for a given magnitude in a given wave band using STMAG definition

    Parameters:
        mag (float): magnitude in AB system
        wvl (float): central wavelength in nanometer
        wvlBand (float): length of wavelength band in nanometres
    Returns:
        float: flux of photons (ph/m2/s)
    """
    return 10 ** (-mag / 2.5) * 10 ** (-21.10/2.5 - 7 + 4 + 8) * c.h.value * c.c.value * wvlBand * wvl


def get_magV(star_type, mag, band):
    """
    Translates magnitude in the specified band in V band

    Parameters:
        star_type (str): only the first letter and number
        mag (float): magnitude in the band defined by "band"
        band (str): band in which the magnitude is expressed, can be B, J, H, K, L, M, N, R, I

    Returns:
        float: magnitude in V band
    """
    if band == 'B':
        magV = mag - colors_dict[star_type[:2] + '.0'][1]
    else:
        ind_col = numpy.where(numpy.asarray(colors_columns) == '(V-' + band + ')')[0][0]
        magV = colors_dict[star_type[:2] + '.0'][ind_col] + mag
    return magV


def flux_star_type(mag, band, star_type, wvl, wvlBand):
    """
    Calculates the flux over a custom wave band, for a given star of a given magnitude in one of the standard bands (I, H, J, etc.)

    Parameters:
        mag (float): magnitude in the band defined by "band"
        band (str): band in which the magnitude is expressed, can be B, J, H, K, L, M, N, R, I
        star_type (str): must be the full star type, i.e. "A0V"
        wvl (float): central wavelength in microns
        wvlBand (float): length of wavelength band in microns

    Returns:
        float: flux of photons (ph/m2/s)
    """
    # get star spectrum, will correspond to a V=0
    Vizier.ROW_LIMIT = -1
    catalog = Vizier.get_catalogs('J/PASP/110/863/synphot')[0]
    pickles_index = numpy.where(catalog["SpType"] == star_type)[0][0] + 1

    ftp = FTP('ftp.stsci.edu')
    ftp.login()

    ss = BytesIO()
    message = ftp.retrbinary('RETR cdbs/grid/pickles/dat_uvk/pickles_uk_' + str(pickles_index) +
                             '.ascii', ss.write)
    ftp.quit()

    # spectrum is in erg/s/cm2/A, first column are wavelengths in A, second column is the spectrum
    erg_spec = numpy.asarray([numpy.asarray(row.split('  '), dtype='float')
                              for row in str(ss.getvalue()).split("\\n")[39:-1]], dtype='float')
    wavelengths = erg_spec[:, 0] * 1e-4

    # translate spectrum in photon/s/m2/micron
    phot_spec = erg_spec[:, 1] * wavelengths / (c.h.value * c.c.value) * 1e-5

    # translate magnitude of band in V
    magV = get_magV(star_type, mag, band)

    # normalise spectrum to correspond to that magnitude
    fac = 10 ** (-magV / 2.5)
    norm_spec = phot_spec * fac

    # integrate over band
    ind = numpy.where((wavelengths >= wvl - wvlBand/2) & (wavelengths <= wvl + wvlBand/2))
    selection = norm_spec[ind]
    flux = (numpy.sum(selection) - (selection[0] + selection[-1]) / 2) * \
           numpy.mean(numpy.diff(wavelengths[ind]))
    return flux


def photons_per_band(mag, mask, pxlScale, expTime, waveband='V'):
        '''
        Calculates the photon flux for a given aperture, star magnitude and wavelength band

        Parameters:
            mag (float): Star apparent magnitude
            mask (ndarray): 2-d pupil mask array, 1 is transparent, 0 opaque
            pxlScale (float): size in metres of each pixel in mask
            expTime (float): Exposure time in seconds
            waveband (string): Waveband

        Returns:
            float: number of photons
        '''

        #Area defined m
        area = mask.sum() * pxlScale**2

        # Flux density photons s^-1 m^-2
        flux_photons = magnitude_to_flux(mag,waveband)

        # Total photons
        photons = flux_photons * expTime * area

        photons = float(photons)

        return photons


def magnitude_to_flux(magnitude, waveband='V'):
    """
    Converts apparent magnitude to a flux of photons
    
    Parameters:
        magnitude (float): Star apparent magnitude
        waveband (string): Waveband of the stellar magnitude, can be U, B, V, R, I, J, H, K, g, r, i, z

    Returns:
        float: Number of photons emitted by the object per second per meter squared

    """

    flux_Jy = FLUX_DICTIONARY[waveband][2] * 10 ** (-0.4 * magnitude)
    flux_photons = flux_Jy * 1.51E7 * FLUX_DICTIONARY[waveband][1] / FLUX_DICTIONARY[waveband][0]  # photons sec^-1 m^-2
    return flux_photons


def flux_to_magnitude(flux, waveband='V'):
    """
    Converts incident flux of photons to the apparent magnitude
    
    Parameters:
        flux (float): Number of photons received from an object per second per meter squared
        waveband (string): Waveband of the measured flux, can be U, B, V, R, I, J, H, K, g, r, i, z

    Returns:
        float: Apparent magnitude
    """

    flux_Jy = flux / (1.51E7 * FLUX_DICTIONARY[waveband][1] / FLUX_DICTIONARY[waveband][0])
    magnitude = float(-2.5 * numpy.log10(flux_Jy / FLUX_DICTIONARY[waveband][2]))
    return magnitude