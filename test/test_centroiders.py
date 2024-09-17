from aotools import image_processing
import numpy
# from nose.tools import raises


def test_centre_of_gravity_single():
    img = numpy.random.random((10, 10))
    com = image_processing.centre_of_gravity(img, 0.1)
    assert(com.shape[0]) == 2


def test_centre_of_gravity_many():
    img = numpy.random.random((5, 10, 10))
    com = image_processing.centre_of_gravity(img, 0.1)
    assert(com.shape[0] == 2)
    assert(com.shape[1] == 5)

def test_centre_of_gravity_value():
    img = numpy.zeros((1, 5, 5))
    img[0, 1:3, 2:4] = 1.
    centroid = image_processing.centre_of_gravity(img)
    numpy.testing.assert_almost_equal(centroid, numpy.array([[2.5], [1.5]]))

def test_brightest_pixel_single():
    img = numpy.random.random((10, 10))
    com = image_processing.brightest_pixel(img, 0.3)
    assert(com.shape[0] == 2)


def test_brightestPxl_many():
    img = numpy.random.random((5, 10, 10))
    com = image_processing.brightest_pixel(img, 0.1)
    assert(com.shape[0] == 2)
    assert(com.shape[1] == 5)


def test_quadCell_single():
    img = numpy.random.random((2, 2))
    com = image_processing.quadCell(img)
    assert(com.shape[0] == 2)


def test_quadCell_many():
    img = numpy.random.random((5, 2, 2))
    com = image_processing.quadCell(img)
    assert(com.shape[0] == 2)
    assert(com.shape[1] == 5)


def test_convolution():
    im = numpy.random.random((10, 10))
    ref = numpy.random.random((10, 10))
    corr = image_processing.cross_correlate(im, ref, padding=1)
    assert(corr.shape == im.shape)


def test_correlation_single():
    im = numpy.random.random((10, 10))
    ref = numpy.random.random((10, 10))
    com = image_processing.correlation_centroid(im, ref, 0.3)
    assert(com.shape[0] == 2)


def test_correlation_many():
    im = numpy.random.random((5, 10, 10))
    ref = numpy.random.random((10, 10))
    com = image_processing.correlation_centroid(im, ref, 0.3)
    assert (com.shape[0] == 2)
    assert(com.shape[1] == 5)


# @raises(ValueError)
# def test_correlation_error():
#     im = numpy.random.random((10))
#     ref = numpy.random.random((10, 10))
#     com = image_processing.correlation_centroid(im, ref, threshold=0.3)



