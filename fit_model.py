"""
This module containts the fitting routine
"""

# General purpose
import numpy as np
import time

#Fitting
import emcee
import pickle

# Astropy
from astropy.io import fits


def make_input_parameters(name,value,minimum,maximum):

    if len(name) != len(minimum) or len(name) != len(maximum) or len(name) != len(value):
        print('Some limits are missing. %d parameteers with %d values, %d minima and %d maxima'
            %(len(name),len(value),len(minimum),len(maximum)))
        return None

    parameters = {}
    for par_name,par_value,par_min,par_max in zip(name,value,minimum,maximum):
        parameters.update({par_name : {'value':par_value,'min':par_min,'max': par_max}})

    return parameters


def update_parameters_values(parameters,value):

    for par_name,par_value in zip(parameters.keys(),value):
        parameters[par_name]['value'] = par_value

    return parameters


def bin_data(data,binning_map):

    if binning_map is not None:

        bin_map = fits.getdata(binning_map)

        if bin_map.shape != data.shape:
            raise Exception('Binning image and observations do not have the same shape (%s, %s)'%(bin_map.shape, data.shape))
    
        else:
            bins = np.unique(bin_map[np.where(bin_map>=0)])
            return np.array([np.nanmean(data[np.where(bin_map == bin_nb)]) for bin_nb in bins])

    else:
        return data

    
def mask_data(obs,mask):

    if mask is not None:

        if np.any(np.unique(mask) != (0,1)):
            raise Exception('Mask should only contain zeros (masked) and ones (valid pixels).')

        elif mask.shape != obs.data.shape:
            raise Exception('Mask and observations do not have the same shape: (%s,%s)'%(mask.shape, Observation.data.shape))

        else:
            mask = mask.astype(float)
            mask[np.where(mask == 0)] = np.nan
            data_obs = obs.data * mask
            data_unc = obs.unc * mask
            return data_obs, data_unc
    else:
        return obs.data,obs.unc



def fit_model(obs,model,parameters,outputname,nsteps=1000,nwalkers=24,mask=None,binning_map=None):

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
    with open(outputname+".pickle",'wb') as f:
        results = {}
        results['chain'] = sampler.chain
        results['lnprobability'] = sampler.lnprobability
        results['parameters'] = parameters
        results['parameters_order'] = parameters.keys()
        results['mask'] = mask
        results['bin_map'] = binning_map
        results['obs'] = obs
        results['model'] = model
        pickle.dump(results,f)

    print('Execution time: %0.4f minutes'%(float(time.time() - start_time)/60))
    
    return results
