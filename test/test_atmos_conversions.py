from aotools.turbulence import atmos_conversions
import numpy


def test_cn2_to_seeing():
    cn2 = 5e-13
    seeing = atmos_conversions.cn2_to_seeing(cn2)
    assert isinstance(seeing, float)


def test_cn2_to_r0():
    cn2 = 5e-13
    r0 = atmos_conversions.cn2_to_r0(cn2)
    assert isinstance(r0, (float, numpy.float64))


def test_r0_to_seeing():
    r0 = 0.1
    seeing = atmos_conversions.r0_to_seeing(r0)
    assert isinstance(seeing, (float, numpy.float64))


def test_conversion_consistency():
    """
    Test internal consistency of conversion by comparing with intermediate steps
    """
    cn2 = 5e-13
    r0 = atmos_conversions.cn2_to_r0(cn2)
    seeing_1 = atmos_conversions.r0_to_seeing(r0)
    seeing_2 = atmos_conversions.cn2_to_seeing(cn2)
    assert seeing_1 == seeing_2


def test_conversion_consistency_inverse():
    """
    Test that the conversions invert correctly
    """
    seeing1 = 1.0
    cn2 = atmos_conversions.seeing_to_cn2(seeing1)
    seeing2 = atmos_conversions.cn2_to_seeing(cn2)
    assert seeing1 == seeing2


def test_coherenceTime():
    cn2 = numpy.array((5e-13, 1e-14))
    v = numpy.array((10, 20))
    tau0 = atmos_conversions.coherenceTime(cn2, v)
    assert isinstance(tau0, (float, numpy.float64))


def test_isoplanaticAngle():
    cn2 = numpy.array((5e-13, 1e-14))
    h = numpy.array((0., 10000.))
    isoplanatic_angle = atmos_conversions.isoplanaticAngle(cn2, h)
    assert isinstance(isoplanatic_angle, (float, numpy.float64))


def test_rytov_variance():
    cn2 = numpy.array((5e-13, 1e-14))
    h = numpy.array((0., 10000.))
    rytov_variance = atmos_conversions.rytov_variance(cn2, h)
    assert isinstance(rytov_variance, (float, numpy.float64))


def test_r0_from_slopes():
    slopes = numpy.random.random((2, 10, 100))
    wavelength = 500e-9
    subapDiam = 0.5
    r0 = atmos_conversions.r0_from_slopes(slopes, wavelength, subapDiam)
    assert isinstance(r0, (float, numpy.float64))


def test_slope_variance_from_r0():
    r0 = 0.15
    wavelength = 500e-9
    subapDiam = 0.5
    variance = atmos_conversions.slope_variance_from_r0(r0, wavelength, subapDiam)
    assert isinstance(variance, (float, numpy.float64))


def test_r0_coversion_consistancy():
    r0 = 0.15
    wavelength = 500e-9
    subapDiam = 0.5
    variance = atmos_conversions.slope_variance_from_r0(r0, wavelength, subapDiam)
    slopes = numpy.random.normal(0., numpy.sqrt(variance), size=(2, 10, 200))
    test_r0 = atmos_conversions.r0_from_slopes(slopes, wavelength, subapDiam)
    numpy.testing.assert_almost_equal(r0, test_r0, 2)
