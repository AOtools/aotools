"""
General equations and definitions describing turbulence statistics
"""
import numpy
from scipy.special import gamma, kv

__all__ = ["phase_covariance","calc_seeing"]

def phase_covariance(r, r0, L0):
    """
    Calculate the phase covariance between two points seperated by `r`, 
    in turbulence with a given `r0 and `L0`.
    Uses equation 5 from Assemat and Wilson, 2006.

    Parameters:
        r (float, ndarray): Seperation between points in metres (can be ndarray)
        r0 (float): Fried parameter of turbulence in metres
        L0 (float): Outer scale of turbulence in metres
    """
    # Make sure everything is a float to avoid nasty surprises in division!
    r = numpy.float32(r)
    r0 = float(r0)
    L0 = float(L0)

    # Get rid of any zeros
    r += 1e-40

    A = (L0 / r0) ** (5. / 3)

    B1 = (2 ** (-5. / 6)) * gamma(11. / 6) / (numpy.pi ** (8. / 3))
    B2 = ((24. / 5) * gamma(6. / 5)) ** (5. / 6)

    C = (((2 * numpy.pi * r) / L0) ** (5. / 6)) * kv(5. / 6, (2 * numpy.pi * r) / L0)

    cov = A * B1 * B2 * C

    return cov


def calc_seeing(r0,lam,l0,r0IsAt500nm=1):
    """Compute seeing from r0, wavelength at which to compute seeing, and L0.
    Note, L0 should be defined at lame.

    Parameters:
       r0 (float, ndarray): Frieds parameter in m.
       lam (float, ndarray): Wavelength in nm.
       l0 (float, ndarray): Outer scale in m.
       r0IsAt500nm (int): Flag, defines where r0 is defined.
    """
    
    lam=lam*1e-9#convert to m
    if type(lam)==type(0.) and lam<1e-10:#probably in m already.
        print "Warning: Wavelength should be defined in nm"
    if r0>2:#probably in cm.  Convert to m.
        print "Warning - r0 might be defined in cm - needs to be in m"

    if r0IsAt500nm:
        r0*=(lam/500e-9)**(6./5)

    seeing = 0.976* lam/r0*180/numpy.pi*3600

    if l0!=0:#outer scale is defined...
        seeing = seeing * numpy.sqrt(1-2.183*(r0/l0)**0.356)
    return seeing
    
