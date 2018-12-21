import numpy as np
import time

import emcee
import pickle

from astropy.io import fits
from .utils import mask_data,bin_data,update_parameters_values


__all__ = ['fit_model']

def fit_model(obs,model,parameters,outputname,nsteps=1000,nwalkers=24,mask=None,binning_map=None):
    """
    Routine that fits the observations using a given model and the emcee sampler.

    We make use of the emcee sampler (http://dfm.io/emcee/current/) to fit the free parameters of
    the model to the observations. We are maximising the following log-probabiluty function:

    $$ln(probability) = ln(priors) + ln(likelihood)$$
        
    with the log likelohood function as:
   
    $$ln(likelihood) = -\frac{1}{2} ( \frac{(data-model)^2}{uncertainty^2} + ln(2 pi uncertainty^2))$$
   
    Both the model and the observations should be instances of the Observations and BaseModel 
    classes from frapy.

    The *parameters* input is a dictionary in the shape:

    parameters = {parameter_name1:{'value':X, 'min':Y, 'max':Z}, parameter_name2:{'value':A, 'min':B, 
    'max':C},...}

    where the parameter_name variables should correspond to the parameters in the model being used;
    value is the starting value of each parameter; and min and max the minimum and maximum values 
    allowed. We assume uniform priors to all parameters (i.e. all values between min and max have
    the same prior probability). Parameters not present in this dictionary will not be varied and 
    will be kept to the value of the input model.

    It is possible to mask part of the data out by using a mask. This should be a 2D array, of the
    same shape as the data, with only zeros (masked values) and ones (valid values). The maximisation
    will be made using only the valid values.

    If the data was binned when deriving the quantity being fit, i.e. if pixels were grouped and
    analysed as a single pixel and that value taken as the value of all the pixels grouped, it is
    possible to include this information using a *binning_map*. This should be a 2D array in which
    pixels of the same bin are given the same value. Pixels in the image that were not analysed (not
    included in the binning) should be given negative values. These are not included in the minimisation.
 
    Parameters
    ----------
    obs: Observation
        An instance of the Observation class
    model: Metallicity_Gradient,Velocity
        A frapy model (based in the BaseModel class)
    parameters: dictionary
        A dictionary containing the parameters of the model to be varied and their limits. Parameters not 
        in this dictionary will not be varied. 
    outputname: str
        Name of the output pickle file.
    nsteps: int
        number of steps of the emcee walkers. Default: 1000
    nwalkers: int
        Number of emcee walkers. Default: 24
    mask: array int
        Array of the same shape as the data containing only zeros (masked values) or ones (valid values). 
        Optional.
    binning_map: array int
        Array of the same shape as the data containing encoding the pixels that were groupped togther.
        Optional.

    Returns
    -------
        Returns a dictionary with:
            sampler chain
            sampler lnprobability
            parameter names in the correct order
            input parameters
            the mask used
            the binning map used 
            the observations 
            the model
        
        This is also saved as a pickles file.
    """

    start_time = time.time()
    
    # Mask data if mask is available
    data_obs, data_unc = mask_data(obs,mask)
   
    # Bin data if voronoi map is available
    data_obs = bin_data(obs.data,binning_map)
    data_unc = bin_data(obs.unc,binning_map)

    # Priors: uniform
    def lnprior(par):
        for par_value,par_name in zip(par,parameters.keys()):
            if par_value < parameters[par_name]['min'] or par_value > parameters[par_name]['max']:
                return -np.inf
        return 0
    

    # Log likelihood function
    def lnprob(current_position,parameters):
        
        lp = lnprior(current_position)
        if lp == -np.inf: 
            return lp
        
        else:
        
            # Create new model
            parameters = update_parameters_values(parameters,current_position)
            model.update_model_parameters(parameters)
            dummy_model = model.make_model()
            convolved_model = model.convolve_with_seeing(seeing=obs.seeing/2.355)

            # Bin Model
            convolved_model = bin_data(convolved_model,binning_map)

            # Calculate likelihood
            inv_sigma2 = 1.0/(data_unc**2)
            lnp = -0.5*(np.nansum((data_obs - convolved_model)**2*inv_sigma2 - np.log(inv_sigma2*np.sqrt(2*np.pi))))
       
            return lnp + lp 
    
    # Prepare initial positions of walkers
    nwalkers = nwalkers
    ndim = len(parameters.keys())
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,args=[parameters])
    initial_position = [ [np.random.uniform(parameters[name]['min'],parameters[name]['max']) for name in parameters.keys() ] for i in range(nwalkers)]
    print('Using %d walkers and fitting %s:'%(nwalkers,parameters.keys()))

    # Fitting
    print('MCMCing for %d steps'%nsteps)
    for i, result in enumerate(sampler.sample(initial_position, iterations=nsteps)):
        if float(i)/nsteps % 0.1 == 0:
            print("%d %%"%(float(i)/nsteps * 100.))
            
    # Save results
    results = {}
    results['chain'] = sampler.chain
    results['lnprobability'] = sampler.lnprobability
    results['parameters_order'] = list(parameters.keys())
    results['mask'] = mask
    results['bin_map'] = binning_map
    results['obs'] = obs
    results['model'] = model       
    with open(outputname+".pickle",'wb') as f:
        pickle.dump(results,f)

    print('Execution time: %0.4f minutes'%(float(time.time() - start_time)/60))
    
    return results
