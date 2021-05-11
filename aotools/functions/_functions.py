import numpy
from numpy import ndarray
from typing import Optional, Tuple, Union


def gaussian2d(size: Union[Tuple[int, int], int], width: Union[Tuple[int, int], int], amplitude: float=1., cent: Optional[Union[ndarray, Tuple[int, int]]]=None) -> ndarray:
    '''
    Generates 2D gaussian distribution


    Args:
        size (tuple, float): Dimensions of Array to place gaussian (y, x)
        width (tuple, float): Width of distribution.
                                Accepts tuple for x and y values in order (y, x).
        amplitude (float, optional): Amplitude of guassian distribution. default is 1.
        cent (tuple, optional): Centre of distribution on grid in order (y, x). Default is middle
    '''

    try:
        ySize = size[0]
        xSize = size[1]
    except (TypeError, IndexError):
        xSize = ySize = size

    try:
        yWidth = float(width[0])
        xWidth = float(width[1])
    except (TypeError, IndexError):
        xWidth = yWidth = float(width)
    
    # If a centre point not given, centre is centre of array
    if cent is None:
        xCent = xSize/2.
        yCent = ySize/2.
    else:
        yCent = cent[0]
        xCent = cent[1]

    X, Y = numpy.meshgrid(range(0, xSize), range(0, ySize))

    image = amplitude * numpy.exp(
        -(((xCent - X) / xWidth) ** 2 + ((yCent - Y) / yWidth) ** 2) / 2)

    return image