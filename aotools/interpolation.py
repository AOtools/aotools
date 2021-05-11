import numpy
from scipy.interpolate import interp2d,RectBivariateSpline
from numpy import ndarray
from typing import Tuple

#a lookup dict for interp2d order (expressed as 'kind')
INTERP_KIND = {1: 'linear', 3:'cubic', 5:'quintic'}


def zoom(array: ndarray, newSize: Tuple[int, int], order: int=3) -> ndarray:
    """
    A Class to zoom 2-dimensional arrays using interpolation

    Uses the scipy `Interp2d` interpolation routine to zoom into an array. Can cope with real of complex data.

    Parameters:
        array (ndarray): 2-dimensional array to zoom
        newSize (tuple): the new size of the required array
        order (int, optional): Order of interpolation to use (1, 3, 5). default is 3

    Returns:
        ndarray : zoom array of new size.
    """
   
    if order not in INTERP_KIND:
       raise ValueError("Order can either be 1, 3, or 5 only")

    try:
        xSize = newSize[0]
        ySize = newSize[1]

    except (IndexError, TypeError):
        xSize = ySize = newSize

    coordsX = numpy.linspace(0, array.shape[0]-1, xSize)
    coordsY = numpy.linspace(0, array.shape[1]-1, ySize)

    #If array is complex must do 2 interpolations
    if array.dtype==numpy.complex64 or array.dtype==numpy.complex128:

        realInterpObj = interp2d(   numpy.arange(array.shape[0]),
                numpy.arange(array.shape[1]), array.real, copy=False, 
                kind=INTERP_KIND[order])
        imagInterpObj = interp2d(   numpy.arange(array.shape[0]),
                numpy.arange(array.shape[1]), array.imag, copy=False,
                kind=INTERP_KIND[order])                 
        return (realInterpObj(coordsY,coordsX) 
                            + 1j*imagInterpObj(coordsY,coordsX))

        

    else:

        interpObj = interp2d(   numpy.arange(array.shape[0]),
                numpy.arange(array.shape[1]), array, copy=False,
                kind=INTERP_KIND[order])

        #return numpy.flipud(numpy.rot90(interpObj(coordsY,coordsX)))
        return interpObj(coordsY,coordsX) 

def zoom_rbs(array: ndarray, newSize: Tuple[int, int], order: int=3) -> ndarray:
    """
    A Class to zoom 2-dimensional arrays using RectBivariateSpline interpolation

    Uses the scipy ``RectBivariateSpline`` interpolation routine to zoom into an array. Can cope with real of complex data. May be slower than above ``zoom``, as RBS routine copies data.

    Parameters:
        array (ndarray): 2-dimensional array to zoom
        newSize (tuple): the new size of the required array
        order (int, optional): Order of interpolation to use. default is 3

    Returns:
        ndarray : zoom array of new size.
    """
    try:
        xSize = newSize[0]
        ySize = newSize[1]

    except IndexError:
        xSize = ySize = newSize

    coordsX = numpy.linspace(0, array.shape[0]-1, xSize)
    coordsY = numpy.linspace(0, array.shape[1]-1, ySize)

    #If array is complex must do 2 interpolations
    if array.dtype==numpy.complex64 or array.dtype==numpy.complex128:
        realInterpObj = RectBivariateSpline(   
                numpy.arange(array.shape[0]), numpy.arange(array.shape[1]), 
                array.real, kx=order, ky=order)
        imagInterpObj = RectBivariateSpline(   
                numpy.arange(array.shape[0]), numpy.arange(array.shape[1]), 
                array.imag, kx=order, ky=order)
                         
        return (realInterpObj(coordsY,coordsX)
                            + 1j*imagInterpObj(coordsY,coordsX))
            
    else:

        interpObj = RectBivariateSpline(   numpy.arange(array.shape[0]),
                numpy.arange(array.shape[1]), array, kx=order, ky=order)


        return interpObj(coordsY,coordsX)
        
        
def binImgs(data: ndarray, n: int) -> ndarray:
    '''
    Bins one or more images down by the given factor
    bins. n must be a factor of data.shape, who knows what happens
    otherwise......
    '''
    shape = numpy.array( data.shape )
    
    n = int(numpy.round(n))

    if len(data.shape)==2:
        shape[-1]/=n
        binnedImgTmp = numpy.zeros( shape, dtype=data.dtype )
        for i in range(n):
            binnedImgTmp += data[:,i::n]
        shape[-2]/=n
        binnedImg = numpy.zeros( shape, dtype=data.dtype )
        for i in range(n):
            binnedImg += binnedImgTmp[i::n,:]

        return binnedImg
    else:
        shape[-1]/=n
        binnedImgTmp = numpy.zeros ( shape, dtype=data.dtype )
        for i in range(n):
            binnedImgTmp += data[...,i::n]

        shape[-2] /= n
        binnedImg = numpy.zeros( shape, dtype=data.dtype )
        for i in range(n):
            binnedImg += binnedImgTmp[...,i::n,:]

        return binnedImg