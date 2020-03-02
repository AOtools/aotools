from aotools import astronomy, circle


def test_photons_per_mag():
    mask = circle(2, 5)
    photons = astronomy.photons_per_mag(5.56, mask, 0.5, 0.3, 10)
    assert type(photons) == float


def test_flux_to_magnitude():
    magnitude = astronomy.flux_to_magnitude(52504716., 'V')
    assert type(magnitude) == float


def test_magnitude_to_flux():
    flux = astronomy.magnitude_to_flux(5.56, 'V')
    assert type(flux) == float


def test_magnitude_to_flux_and_flux_to_magnitude():
    flux = astronomy.magnitude_to_flux(5.56, 'V')
    magnitude = astronomy.flux_to_magnitude(flux, 'V')
    assert magnitude == 5.56


def test_photons_per_band():
    photons = astronomy.photons_per_band(5.56, circle(2, 5), 0.5, 0.001)
    assert type(photons) == float


def test_flux_AB_mag():
    flux = astronomy.flux_AB_mag(10, 0.91, 0.13)
    assert type(flux) == float


def test_flux_STMAG_mag():
    flux = astronomy.flux_STMAG_mag(10, 0.91, 0.13)
    assert type(flux) == float


def test_get_magV():
    magV = astronomy.get_magV('A0', 0, 'H')
    assert type(magV) == float
    assert magV == -0.19


def test_flux_star_type():
    flux = astronomy.flux_star_type(10, 'H', 'A0V', 1.1, 0.2)
    assert type(flux) == astronomy.numpy.float64

