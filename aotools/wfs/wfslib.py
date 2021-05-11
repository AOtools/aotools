"""
A library of functions which may be of use to analyse WFS data
"""

import numpy
from numpy import ndarray
from typing import Tuple, Union

# Best range for python2 and 3
try:
    range = xrange
except:
    pass


def findActiveSubaps(subaps: int, mask: ndarray, threshold: float, returnFill: bool=False) -> Union[ndarray, Tuple[ndarray, ndarray]]:
    """
    Finds the subapertures which are "seen" be through the
    pupil function. Returns the coords of those subapertures

    Parameters:
        subaps (int): The number of subaps in x (assumes square)
        mask (ndarray): A pupil mask, where is transparent when 1, and opaque when 0
        threshold (float): The mean value across a subap to make it "active"
        returnFill (optional, bool): Return an array of fill-factors

    Returns:
        ndarray: An array of active subap coords
    """

    subapCoords = []
    xSpacing = mask.shape[0]/float(subaps)
    ySpacing = mask.shape[1]/float(subaps)

    if returnFill:
        fills = []

    for x in range(subaps):
        for y in range(subaps):
            subap = mask[
                    int(numpy.round(x*xSpacing)): int(numpy.round((x+1)*xSpacing)),
                    int(numpy.round(y*ySpacing)): int(numpy.round((y+1)*ySpacing))
                    ]

            if subap.mean() >= threshold:
                subapCoords.append( [x*xSpacing, y*ySpacing])
                if returnFill:
                    fills.append(subap.mean())

    subapCoords = numpy.array( subapCoords )

    if returnFill:
        return subapCoords, numpy.array(fills)
    else:
        return subapCoords


def computeFillFactor(mask: ndarray, subapPos: ndarray, subapSpacing: int) -> ndarray:
    """
    Calculate the fill factor of a set of sub-aperture co-ordinates with a given
    pupil mask.

    Parameters:
        mask (ndarray): Pupil mask
        subapPos (ndarray): Set of n sub-aperture co-ordinates (n, 2)
        subapSpacing: Number of mask pixels between sub-apertures

    Returns:
        list: fill factor of sub-apertures
    """

    fills = numpy.zeros(len(subapPos))
    for i, (x, y) in enumerate(subapPos):
        x1 = int(round(x))
        x2 = int(round(x + subapSpacing))
        y1 = int(round(y))
        y2 = int(round(y + subapSpacing))
        fills[i] = mask[x1:x2, y1:y2].mean()

    return fills

def make_subaps_2d(data: ndarray, mask: ndarray) -> ndarray:
    """
    Fills in a pupil shape with 2-d sub-apertures

    Parameters:
        data (ndarray): slope data, of shape, (frames, 2, nSubaps)
        mask (ndarray): 2-d array of shape (nxSubaps, nxSubaps), where 1 indicates a valid subap position and 0 a masked subap 
    """ 
    n_frames = data.shape[0]
    n_subaps = data.shape[-1]
    nx_subaps = mask.shape[0]

    subaps_2d = numpy.zeros((data.shape[0], 2, mask.shape[0], mask.shape[1]), dtype=data.dtype)

    n_subap = 0
    for x in range(nx_subaps):
        for y in range(nx_subaps):
            if mask[x, y] == 1:
                subaps_2d[:, :, x, y] = data[:, :, n_subap]
                n_subap += 1

    return subaps_2d
