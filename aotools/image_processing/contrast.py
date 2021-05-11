"""
Contrast
--------

Functions for calculating the contrast of an image.
"""
from numpy import ndarray

def image_contrast(image: ndarray) -> float:
    """
    Calculates the 'Michelson' contrast.

    Uses a method by Michelson (Michelson, A. (1927). Studies in Optics. U. of Chicago Press.), to calculate the contrast ratio of an image. Uses the formula:
        (img_max - img_min)/(img_max + img_min)

    Parameters:
        image (ndarray): Image array

    Returns:
        float: Contrast value
    """

    contrast = (image.max() - image.min()) / (image.max() + image.min())

    return float(contrast)


def rms_contrast(image: ndarray) -> float:
    """
    Calculates the RMS contrast - basically the standard deviation of the image

    Parameters:
        image (ndarray): Image array

    Returns:
        float: Contrast value
    """

    image /= image.max()

    return float(image.std())