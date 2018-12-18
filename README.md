# FRApy
## Fitting Resolved Arcs with python

## Intro

Gravitational arcs are galaxies in the backgroud of a cluster of galaxies or a massive galaxy. These massive objects act as a *magifying lens*, increasing the size of the backgroud galaxies in the sky, which allows us to resolve them to smaller spatial scales. However, the magnification is not uniform, and the lensed galaxies appear 'distorted' in the sky, typically in an arc-like shape (hence the name graviational 'arcs').

[logo]: https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 2"
![alt text](https://goo.gl/images/iyYaVS "Lensed Galaxy in the cluster Abell 370 (Credits:NASA, ESA, A. Koekemoer, M. Jauzac, C. Steinhardt, and the BUFFALO team)")


The spatial distortion makes the analysis of these objects more difficult, especially when several images of the same object (called multiple images) are formed. FRApy allows these images with


## What you'll need

	1. Displacement maps. COMPLETE THIS
	2. An analytical model that you are trying to fit. You'll have to assume some hypothesis or model of your data (such as a linear metallicity gradient, an arctangent function to fit velocity, or an exponential surface brightness), that will be fit to the data. If you are interested in seeing the data 'as it is in source plane', i.e. in inverting the lensing equations, use LensTool, or similar codes.

## Install

	???

## Example

	FRApy is a python module. A minimum working example would be something like:

	```python
	from obs_and_model import Observation,Metallicity_Gradient
	from fit_model import fit_model,make_input_parameters
	from check_fit import Output

	# Load Observations
	obs = Observation(z=0.611,
                  data_path='Demo_data/AS1063_map_metallicity.fits',
                  unc_path='Demo_data/AS1063_map_metallicity_unc.fits',
                  seeing = 1.03/0.2)
    # Load Model, in this case a linear metallicity gradient
    model = Metallicity_Gradient(zlens=0.322,
              dplx_path='Demo_data/AS1063_dplx.fits',
              dply_path='Demo_data/AS1063_dply.fits')
    model.create_displacement_maps_for_object(obs)

    # Fit the data
    input_par = make_input_parameters(name    = ('cx', 'cy',  'q', 'pa', 'z_grad', 'z_0'),
                                 value   = (  29,   23,  0.7,   20,    -0.02, 9.0),
                                 minimum = (  28,   22,  0.4,  -20,     -0.1, 8.5),
                                 maximum = (  33,   27,  0.9,   90,      0.0, 9.5))
	out = fit_model(obs,model,input_par,'output',nsteps=2000,nwalkers=24)

	# Inspect the fit
	results = Output('as1063_metallicity')
	results.best_parameters()
	```

## More info

	You can find a more detailled example in XXX and XXX.
	These were galaxies analysied in Patricio et al. 2018 and Patricio et al. 2019.

## Authors

	Main author: Vera Patricio
	Lensing Specialist: Johan Richard

## Licence
 ???

