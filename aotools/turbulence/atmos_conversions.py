"""
Atmospheric Parameter Conversions
---------------------------------

Functions for converting between different atmospheric parameters,

"""

import numpy


def cn2_to_seeing(cn2, lamda=500.E-9):
    """
    Calculates the seeing angle from the integrated Cn2 value

    Parameters:
        cn2 (float): integrated Cn2 value in m^2/3
        lamda : wavelength

    Returns:
        seeing angle in arcseconds
    """
    r0 = cn2_to_r0(cn2,lamda)
    seeing = r0_to_seeing(r0,lamda)
    return seeing


def seeing_to_cn2(seeing, lamda=500.E-9):
    """
    Calculates the integrated Cn2 value from the seeing

    Parameters:
        seeing (float): seeing in arcseconds
        lamda : wavelength

    Returns:
        integrated Cn2 value in m^2/3
    """
    r0 = seeing_to_r0(seeing,lamda)
    cn2 = r0_to_cn2(r0,lamda)
    return cn2


def cn2_to_r0(cn2, lamda=500.E-9):
    """
    Calculates r0 from the integrated Cn2 value

    Parameters:
        cn2 (float): integrated Cn2 value in m^2/3
        lamda : wavelength

    Returns:
        r0 in m
    """
    r0=(0.423*(2*numpy.pi/lamda)**2*cn2)**(-3./5.)
    return r0


def r0_to_cn2(r0, lamda=500.E-9):
    """
    Calculates integrated Cn2 value from r0

    Parameters:
        r0 (float): r0 in cm
        lamda : wavelength

    Returns:
        cn2 (float): integrated Cn2 value in m^2/3
    """
    cn2 = r0**(-5./3.)/(0.423*(2*numpy.pi/lamda)**2)
    return cn2


def r0_to_seeing(r0, lamda=500.E-9):
    """
    Calculates the seeing angle from r0

    Parameters:
        r0 (float): Freid's parameter in cm
        lamda : wavelength

    Returns:
        seeing angle in arcseconds
    """
    return (0.98*lamda/r0)*180.*3600./numpy.pi


def seeing_to_r0(seeing, lamda=500.E-9):
    """
    Calculates r0 from seeing

    Parameters:
        seeing (float): seeing angle in arcseconds
        lamda : wavelength

    Returns:
        r0 (float): Freid's parameter in cm
    """
    return 0.98*lamda/(seeing*numpy.pi/(180.*3600.))


def coherenceTime(cn2, v, lamda=500.E-9, axis=-1):
    """
    Calculates the coherence time from profiles of the Cn2 and wind velocity

    Parameters:
        cn2 (array): Cn2 profile in m^2/3. >1D arrays are supported (i.e. multiple
            profiles), in which case the default assumption is that the turbulent 
            layers are along the final (-1) axis. 
        v (array): profile of wind velocity, same altitude scale as cn2 
        lamda : wavelength
        axis (int, optional): axis over which to integrate (for >1D inputs)

    Returns:
        coherence time in seconds
    """
    Jv = (cn2*(v**(5./3.))).sum(axis)
    tau0 = (Jv**(-3./5.))*0.057*lamda**(6./5.)
    return tau0


def isoplanaticAngle(cn2, h, lamda=500.E-9, axis=-1):
    """
    Calculates the isoplanatic angle from the Cn2 profile

    Parameters:
        cn2 (array): Cn2 profile in m^2/3.  >1D arrays are supported (i.e. multiple
            profiles), in which case the default assumption is that the turbulent 
            layers are along the final (-1) axis. 
        h (Array): Altitude levels of cn2 profile in m
        lamda : wavelength
        axis (int, optional): axis over which to integrate (for >1D inputs)

    Returns:
        isoplanatic angle in arcseconds
    """
    Jh = (cn2*(h**(5./3.))).sum(axis)
    iso = 0.057*lamda**(6./5.)*Jh**(-3./5.)*180.*3600./numpy.pi
    return iso


def rytov_variance(cn2, h, lamda=500.E-9, axis=-1):
    """
    Calculates plane wave Rytov variance from the Cn2 profile

    Parameters:
        cn2 (numpy.ndarray): Cn2 profile in m^2/3.  >1D arrays are supported (i.e. multiple
            profiles), in which case the default assumption is that the turbulent 
            layers are along the final (-1) axis. 
        h (numpy.ndarray): Altitude levels of cn2 profile in m
        lamda (float, optional): wavelength, default 500 nm
        axis (int, optional): axis over which to integrate (for >1D inputs)

    Returns:
        Rytov variance (float)
    """
    k = 2. * numpy.pi / lamda
    return 2.25 * k**(7./6.) * (cn2*h**(5./6.)).sum(axis)


def r0_from_slopes(slopes, wavelength, subapDiam):
    """
    Measures the value of R0 from a set of WFS slopes.

    Uses the equation in Saint Jaques, 1998, PhD Thesis, Appendix A to calculate the value of atmospheric seeing parameter, r0, that would result in the variance of the given slopes.

    Parameters:
        slopes (ndarray): A 3-d set of slopes in radians, of shape (dimension, nSubaps, nFrames)
        wavelength (float): The wavelegnth of the light observed
        subapDiam (float) The diameter of each sub-aperture

    Returns:
        float: An estimate of r0 for that dataset.

    """
    slopeVar = slopes.var(axis=(-1))

    r0 = ((0.162 * (wavelength ** 2) * subapDiam ** (-1. / 3)) / slopeVar) ** (3. / 5)

    r0 = r0.mean()

    return r0


def slope_variance_from_r0(r0, wavelength, subapDiam):
    """
    Uses the equation in Saint Jaques, 1998, PhD Thesis, Appendix A to calculate the slope variance resulting from a
    value of r0.

    Parameters:
        r0 (float): Fried papamerter of turubulence in metres
        wavelength (float): Wavelength of light in metres (where 1e-9 is 1nm)
        subapDiam (float): Diameter of the aperture in metres

    Returns:
        The expected slope variance for a given r0 ValueError

    """

    slope_var = 0.162 * (wavelength ** 2) * r0 ** (-5. / 3) * subapDiam ** (-1. / 3)

    return slope_var
