from aotools.functions import ell_zernike
import matplotlib.pyplot as plt
import numpy as np


def test_ZernikeEllipticalaperture():
    # Define parameters for the ZernikeEllipticalaperture instance
    nterms = 6  # Number of Zernike terms
    npix = 256  # Number of pixels in each dimension
    a = 1.0  # Semi-major axis of the elliptical aperture
    b = 0.5  # Semi-minor axis of the elliptical aperture

    zernike_instance = ell_zernike.ZernikeEllipticalaperture(nterms, npix, a, b)

    assert zernike_instance.ell_aperture_mask.shape == (npix, npix), "Aperture mask shape is incorrect"

    assert np.all(np.isin(zernike_instance.ell_aperture_mask, [0, 1])), "Aperture mask should contain only 0s and 1s"

    assert zernike_instance.E.shape == (nterms, npix, npix), "Zernike modes shape is incorrect"

    assert np.any(zernike_instance.E[0][zernike_instance.GenerateEllipticalAperture() == 1]) != 0, "First Zernike mode should have non-zero values in the aperture"

    expected_number_of_modes = nterms
    assert zernike_instance.E.shape[0] == expected_number_of_modes, f"Expected {expected_number_of_modes} Zernike modes, got {zernike_instance.E.shape[0]}"

    phi = zernike_instance.EllZernikeMap()
    assert phi.shape == (npix, npix), "Output shape is incorrect when no coefficients are provided."

    coeff = np.random.random(nterms)
    phi_with_coeff = zernike_instance.EllZernikeMap(coeff)
    assert phi_with_coeff.shape == (npix, npix), "Output shape is incorrect with provided coefficients."


test_ZernikeEllipticalaperture()
