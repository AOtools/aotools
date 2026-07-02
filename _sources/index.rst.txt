.. aotools documentation master file, created by
   sphinx-quickstart on Thu Mar 10 10:09:36 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the AOtools Documentation!
=====================================

Introduction
++++++++++++
AOtools is an attempt to gather together many common tools, equations and functions for Adaptive Optics.
The idea is to provide an alternative to the current model, where a new AO researcher must write their own library
of tools, many of which implement theory and ideas that have been in existance for over 20 years. AOtools will hopefully
provide a common place to put these tools, which through use and bug fixes, should become reliable and well documented.

Installation
------------
AOtools uses mainly standard library functions, and all attempts have been made to avoid adding unneccessary dependencies.
Currently, the library requires only

- numpy
- scipy
- astropy
- matplotlib

Contents
++++++++

.. toctree::
   :maxdepth: 2

   introduction
   turbulence
   image_processing
   zernike
   pupil
   wfs
   opticalpropagation
   astronomy
   fft
   interpolation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
