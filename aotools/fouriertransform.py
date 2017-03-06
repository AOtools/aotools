"""
Module containing useful FFT based function and classes
"""
import numpy

def ft(data, delta):
    """
    A properly scaled 1-D FFT

    Parameters:
        data (ndarray): An array on which to perform the FFT
        delta (float): Spacing between elements

    Returns:
        ndarray: scaled FFT
    """
    DATA = numpy.fft.fftshift(
            numpy.fft.fft(
                    numpy.fft.fftshift(data, axes=(-1))),
            axes=(-1)) * delta
    return DATA

def ift(DATA, delta_f):
    """
    Scaled inverse 1-D FFT

    Parameters:
        DATA (ndarray): Data in Fourier Space to transform
        delta_f (ndarray): Frequency spacing of grid

    Returns:
        ndarray: Scaled data in real space
    """

    data = numpy.fft.ifftshift(
            numpy.fft.ifft(
                    numpy.fft.ifftshift(DATA, axes=(-1))),
            axes=(-1)) * len(DATA) * delta_f

    return data


def ft2(data, delta):
    """
    A properly scaled 2-D FFT

    Parameters:
        data (ndarray): An array on which to perform the FFT
        delta (float): Spacing between elements

    Returns:
        ndarray: scaled FFT
    """

    DATA = numpy.fft.fftshift(
            numpy.fft.fft2(
                    numpy.fft.fftshift(data, axes=(-1,-2))
                    ), axes=(-1,-2)
            )*delta**2

    return DATA

def ift2(DATA, delta_f):
    """
    Scaled inverse 2-D FFT

    Parameters:
        DATA (ndarray): Data in Fourier Space to transform
        delta_f (ndarray): Frequency spacing of grid

    Returns:
        ndarray: Scaled data in real space
    """
    N = DATA.shape[0]
    g = numpy.fft.ifftshift(
            numpy.fft.ifft2(
                    numpy.fft.ifftshift(DATA))) * (N * delta_f)**2

    return g

def rft(data, delta):
    """
    A properly scaled real 1-D FFT

    Parameters:
        data (ndarray): An array on which to perform the FFT
        delta (float): Spacing between elements

    Returns:
        ndarray: scaled FFT
    """
    DATA = numpy.fft.fftshift(
            numpy.fft.rfft(
                    numpy.fft.fftshift(data, axes=(-1))),
            axes=(-1)) * delta
    return DATA

def irft(DATA, delta_f):
    """
    Scaled real inverse 1-D FFT

    Parameters:
        DATA (ndarray): Data in Fourier Space to transform
        delta_f (ndarray): Frequency spacing of grid

    Returns:
        ndarray: Scaled data in real space
    """

    data = numpy.fft.ifftshift(
            numpy.fft.irfft(
                    numpy.fft.ifftshift(DATA, axes=(-1))),
            axes=(-1)) * len(DATA) * delta_f

    return data

def rft2(data, delta):
    """
    A properly scaled, real 2-D FFT

    Parameters:
        data (ndarray): An array on which to perform the FFT
        delta (float): Spacing between elements

    Returns:
        ndarray: scaled FFT
    """
    DATA = numpy.fft.fftshift(
            numpy.fft.rfft2(
                    numpy.fft.fftshift(data, axes=(-1,-2))
                    ), axes=(-1,-2)
            )*delta**2

    return DATA

def irft2(DATA, delta_f):
    """
    Scaled inverse real 2-D FFT

    Parameters:
        DATA (ndarray): Data in Fourier Space to transform
        delta_f (ndarray): Frequency spacing of grid

    Returns:
        ndarray: Scaled data in real space
    """
    data = numpy.fft.ifftshift(
            numpy.fft.irfft2(
                    numpy.fft.ifft2(DATA), s=(-1,-2)
                    )
            ) * (DATA.shape[0]*delta_f)**2
    return data
