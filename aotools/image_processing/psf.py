"""
PSF Analysis
------------

Functions for analysing PSFs.

"""

import numpy
from .. import functions
from numpy import ndarray
from typing import Optional, Tuple, Union


def azimuthal_average(data: ndarray) -> ndarray:
    """
    Measure the azimuthal average of a 2d array

    Args:
        data (ndarray): A 2-d array of data

    Returns:
        ndarray: A 1-d vector of the azimuthal average
    """
    size = data.shape[0]
    avg = numpy.empty(int(size / 2), dtype="float")
    for i in range(int(size / 2)):
        ring = functions.pupil.circle(i + 1, size) - functions.pupil.circle(i, size)
        avg[i] = (ring * data).sum() / (ring.sum())

    return avg


def encircled_energy(data: ndarray,
                    fraction: float=0.5, center: None=None,
                    eeDiameter: bool=True) -> Union[float, Tuple[ndarray, ndarray]]:
    """
        Return the encircled energy diameter for a given fraction
        (default is ee50d).
        Can also return the encircled energy function.
        Translated and extended from YAO.

        Parameters:
            data : 2-d array
            fraction : energy fraction for diameter calculation
                default = 0.5
            center : default = center of image
            eeDiameter : toggle option for return.
                If False returns two vectors: (x, ee(x))
                Default = True
        Returns:
            Encircled energy diameter
            or
            2 vectors: diameters and encircled energies

    """
    dim = data.shape[0] // 2
    if center is None:
        center = [dim, dim]
    xc = center[0]
    yc = center[1]
    e = 1.9
    npt = 20
    rad = numpy.linspace(0, dim**(1. / e), npt)**e
    ee = numpy.empty(rad.shape)

    for i in range(npt):
        pup = functions.pupil.circle(rad[i],
                           int(dim) * 2,
                           circle_centre=(xc, yc),
                           origin='corner')
        rad[i] = numpy.sqrt(numpy.sum(pup) * 4 / numpy.pi)  # diameter
        ee[i] = numpy.sum(pup * data)

    rad = numpy.append(0, rad)
    ee = numpy.append(0, ee)
    ee /= numpy.sum(data)
    xi = numpy.linspace(0, dim, int(4 * dim))
    yi = numpy.interp(xi, rad, ee)

    if eeDiameter is False:
        return xi, yi
    else:
        ee50d = float(xi[numpy.argmin(numpy.abs(yi - fraction))])
        return ee50d