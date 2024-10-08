from aotools import functions
import numpy
import matplotlib.pyplot as plt


def test_zernIndex():
    results = ([0, 0], [1, 1], [1, -1], [2, 0], [2, -2], [2, 2])
    for i in range(1, 6):
        index = functions.zernIndex(i)
        assert(index == results[i-1])

def test_makegammas():
    gammas = functions.makegammas(5)
    assert(gammas.shape == (2, 21, 21))


def test_zenikeRadialFunc():
    coordinates = numpy.linspace(-1, 1, 32)
    x, y = numpy.meshgrid(coordinates, coordinates)
    r = numpy.sqrt(x ** 2 + y ** 2)

    radial_function = functions.zernikeRadialFunc(5, 3, r)
    # test no casting error for high order modes
    _ = functions.zernikeRadialFunc(21, 1, r)
    _ = functions.zernikeRadialFunc(50, 16, r)
    assert(radial_function.shape == (32, 32))


def test_zernike_nm():
    zernike_array = functions.zernike_nm(6, 2, 32)
    assert(zernike_array.shape == (32, 32))


def test_zernike():
    zernike_array = functions.zernike_noll(9, 32)
    assert(zernike_array.shape == (32, 32))


def test_zernikeArray_single():
    zernike_array = functions.zernikeArray(10, 32)
    assert(zernike_array.shape == (10, 32, 32))
    plt.imshow(zernike_array[0])
    plt.show()


def test_zernikeArray_list():
    zernike_array = functions.zernikeArray([2, 3, 4], 32)
    assert(zernike_array.shape == (3, 32, 32))


def test_zernikeArray_comparison():
    full_zernike_array = functions.zernikeArray(10, 32)
    subset_zernike_array = functions.zernikeArray([2, 3, 4], 32)
    assert(numpy.allclose(full_zernike_array[1:4], subset_zernike_array))


def test_phaseFromZernikes():
    phase_map = functions.phaseFromZernikes([1, 2, 3, 4, 5], 32)
    assert(phase_map.shape == (32, 32))


