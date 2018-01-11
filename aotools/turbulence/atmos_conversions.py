import numpy

def cn2_to_seeing(cn2,lamda=500.E-9):
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

def cn2_to_r0(cn2,lamda=500.E-9):
    """
    Calculates r0 from the integrated Cn2 value

    Parameters:
        cn2 (float): integrated Cn2 value in m^2/3
        lamda : wavelength

    Returns:
        r0 in cm
    """
    r0=(0.423*(2*numpy.pi/lamda)**2*cn2)**(-3./5.)
    return r0

def r0_to_seeing(r0,lamda=500.E-9,L0=None,r0IsAt500nm=True):
    """
    Calculates the seeing angle from r0 and L0 (optionally)

    Parameters:
        r0 (float): Freid's parameter in m
        lamda : wavelength in m
        L0 (float): Outer scale in m.
        r0IsAt500nm (Flag): Whether r0 is defined at 500nm.

    Returns:
        seeing angle in arcseconds
    """

    
    if type(lam)==type(0.) and lam>1:#probably in nm
        print( "Warning: Wavelength should be defined in m")
    if r0>2:#probably in cm.  Convert to m.
        print( "Warning - r0 might be defined in cm - needs to be in m")

    if r0IsAt500nm:
        r0*=(lam/500e-9)**(6./5)

    seeing = 0.976* lam/r0*180/numpy.pi*3600

    if l0!=0:#outer scale is defined...
        seeing = seeing * numpy.sqrt(1-2.183*(r0/l0)**0.356)
    return seeing
    

    #return (0.98*lamda/r0)*180.*3600./numpy.pi

def coherenceTime(cn2,v,lamda=500.E-9):
    """
    Calculates the coherence time from profiles of the Cn2 and wind velocity

    Parameters:
        cn2 (array): Cn2 profile in m^2/3
        v (array): profile of wind velocity, same altitude scale as cn2 
        lamda : wavelength

    Returns:
        coherence time in seconds
    """
    Jv = (cn2*(v**(5./3.))).sum()
    tau0 = float((Jv**(-3./5.))*0.057*lamda**(6./5.))
    return tau0

def isoplanaticAngle(cn2,h,lamda=500.E-9):
    """
    Calculates the isoplanatic angle from the Cn2 profile

    Parameters:
        cn2 (array): Cn2 profile in m^2/3
        h (Array): Altitude levels of cn2 profile in m
        lamda : wavelength

    Returns:
        isoplanatic angle in arcseconds
    """
    Jh = (cn2*(h**(5./3.))).sum()
    iso = float(0.057*lamda**(6./5.)*Jh**(-3./5.)  *180.*3600./numpy.pi)
    return iso
