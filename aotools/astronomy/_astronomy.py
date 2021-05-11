import numpy
from numpy import ndarray

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


def photons_per_mag(mag: float, mask: ndarray, pixel_scale: float, wvlBand: float, exposure_time: int) -> float:
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


def photons_per_band(mag: float, mask: ndarray, pxlScale: float, expTime: float, waveband: str='V') -> float:
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


def magnitude_to_flux(magnitude: float, waveband: str='V') -> float:
    """
    Converts apparent magnitude to a flux of photons
    
    Parameters:
        magnitude (float): Star apparent magnitude
        waveband (string): Waveband of the stellar magnitude, can be U, B, V, R, I, J, H, K, g, r, i, z

    Returns:
        float: Number of photons emitted by the object per second per meter squared

    """

    flux_Jy = FLUX_DICTIONARY[waveband][2] * 10 ** (-0.4 * magnitude)
    flux_photons = flux_Jy * 1.51E7 * FLUX_DICTIONARY[waveband][1]  # photons sec^-1 m^-2
    return flux_photons


def flux_to_magnitude(flux: float, waveband: str='V') -> float:
    """
    Converts incident flux of photons to the apparent magnitude
    
    Parameters:
        flux (float): Number of photons received from an object per second per meter squared
        waveband (string): Waveband of the measured flux, can be U, B, V, R, I, J, H, K, g, r, i, z

    Returns:
        float: Apparent magnitude
    """

    flux_Jy = flux / (1.51E7 * FLUX_DICTIONARY[waveband][1])
    magnitude = float(-2.5 * numpy.log10(flux_Jy / FLUX_DICTIONARY[waveband][2]))
    return magnitude