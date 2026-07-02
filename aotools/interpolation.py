import numpy
from scipy.interpolate import RectBivariateSpline
#a lookup dict for interp2d order (expressed as 'kind')
INTERP_KIND = {1: 'linear', 3:'cubic', 5:'quintic'}


def zoom_rbs(array, newSize, order=3):
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
        
# for compatibility with new Scipy where interp2d is dead...
zoom = zoom_rbs

def binImgs(data, n):
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