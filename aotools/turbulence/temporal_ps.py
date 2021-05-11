"""
Temporal Power Spectra
----------------------

Turbulence gradient temporal power spectra calculation and plotting

:author: Andrew Reeves
:date: September 2016
"""

import numpy
from matplotlib import pyplot
from scipy import optimize
from numpy import ndarray
from typing import Tuple


def calc_slope_temporalps(slope_data: ndarray) -> Tuple[ndarray, ndarray]:
    """
    Calculates the temporal power spectra of the loaded centroid data.

    Calculates the Fourier transform over the number frames of the centroid
    data, then returns the square of the  mean of all sub-apertures, for x
    and y. This is a temporal power spectra of the slopes, and should adhere
    to a -11/3 power law for the slopes in the wind direction, and -14/3 in
    the direction tranverse to the wind direction. See Conan, 1995 for more.

    The phase slope data should be split into X and Y components, with all X data first, then Y (or visa-versa).

    Parameters:
        slope_data (ndarray): 2-d array of shape (..., nFrames, nCentroids)

    Returns:
        ndarray: The temporal power spectra of the data for X and Y components.
    """

    n_frames = slope_data.shape[-2]

    # Only take half result, as FFT mirrors
    tps = abs(numpy.fft.fft(slope_data, axis=-2)[..., :int(n_frames/2), :])**2

    # Find mean across all sub-aps
    tps = (abs(tps)**2)
    mean_tps = tps.mean(-1)
    tps_err = tps.std(-1)/numpy.sqrt(tps.shape[-1])

    return mean_tps, tps_err


def get_tps_time_axis(frame_rate: int, n_frames: int) -> ndarray:
    """
    Parameters:
        frame_rate (float): Frame rate of detector observing slope gradients (Hz)
        n_frames: (int): Number of frames used for temporal power spectrum

    Returns:
        ndarray: Time values for temporal power spectra plots
    """

    t_vals = numpy.fft.fftfreq(n_frames, 1./frame_rate)[:int(n_frames/2)]

    return t_vals


def plot_tps(slope_data, frame_rate):
    """
    Generates a plot of the temporal power spectrum/a for a data set of phase gradients

    Parameters:
        slope_data (ndarray):  2-d array of shape (..., nFrames, nCentroids)
        frame_rate (float): Frame rate of detector observing slope gradients (Hz)

    Returns:
        tuple: The computed temporal power spectrum/a, and the time axis data

    """
    n_frames = slope_data.shape[-2]

    tps, tps_err = calc_slope_temporalps(slope_data)

    t_axis_data = get_tps_time_axis(frame_rate, n_frames)

    fig = pyplot.figure()
    ax = fig.add_subplot(111)

    # plot each power spectrum
    for i_ps, ps in enumerate(tps):
        ax.loglog(t_axis_data, ps, label="Spectrum {}".format(i_ps))

    ax.legend()

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power (arbitrary units)")

    pyplot.show()

    return tps, tps_err, t_axis_data

def fit_tps(tps, t_axis, D, V_guess=10, f_noise_guess=20, A_guess=9, tps_err=None, plot=False):
    """
    Runs minimization routines to get t0.

    Parameters:
        tps (ndarray): The temporal power spectrum to fit
        axis (str): fit parallel ('par') or perpendicular ('per') to wind direction
        D (float): (Sub-)Aperture diameter

    Returns:
    """
    if plot:
        fig = pyplot.figure()

    opt_result = optimize.minimize(
            test_tps_fit_minimize_func,
            x0=(V_guess, f_noise_guess, A_guess),
            args=(tps, t_axis, D, tps_err, plot),
            method="COBYLA")

    print(opt_result)
