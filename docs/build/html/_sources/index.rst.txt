
Welcome to FRApy's documentation!
=================================

FRApy - Fitting Resolved Arcs with python - is python code designed to fit gravitational arcs in image plane with analytical modes, such as metallicity gradients and velocity fields, taking into account the lensing distortions.


Intro
=================================

Graviatational arcs
**************************

Gravitational arcs are galaxies in the backgroud of a cluster of galaxies or a massive galaxy. These massive objects act as a *magifying lens*, increasing the size of the backgroud galaxies in the sky, which allows us to resolve them to smaller spatial scales. However, the magnification is not uniform, and the lensed galaxies appear 'distorted' in the sky, typically in an arc-like shape (hence the name graviational 'arcs').

PUT IMAGE

The spatial distortion makes the analysis of these objects more difficult, especially when several images of the same object (called multiple images) are formed. FRApy allows these images with


What you'll need
**************************

1. Displacement maps. COMPLETE THIS

2. An analytical model that you are trying to fit. You'll have to assume some hypothesis or model of your data (such as a linear metallicity gradient, an arctangent function to fit velocity, or an exponential surface brightness), that will be fit to the data. If you are interested in seeing the data 'as it is in source plane', i.e. in inverting the lensing equations, use LensTool, or similar codes.



Install
=================================

Add dependencies

Documentation 
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Observations
**************************
.. automodule:: frapy.observations
	:members:

Models
**************************
.. automodule:: frapy.models
	:members:

Fitting
**************************
.. automodule:: frapy.fit_model
	:members:

Explore the output
**************************
.. automodule:: frapy.check_fit
	:members:

Miscelaneous
**************************
.. automodule:: frapy.utils
	:members:




Authors and References
======================



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
