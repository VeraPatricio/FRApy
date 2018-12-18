"""
This module containts miscelaneous functions used in the fitting.
"""

import numpy as np

from astropy.io import fits

__all__ = ['make_input_parameters']


def make_input_parameters(name,value,minimum,maximum):
    """ 
    Outputs a parameter dictionary to be used in the fit_model function.

    This dictionary has the following form:

        {parameter_name1:{'value':X, 'min':Y, 'max':Z},
        parameter_name2:{'value':A, 'min':B, 'max':C},
        ...
        }

    Parameters
    ----------
    name: array str
        An array of strings containing the names of the model's parameters to be fitted.

    value: array float
        The initial values of these parameters (ACTUALLY NOT USED)
    minimum: array float
        The minimum value allowed for each parameter
    maximum: array float
        The maximum value allowed for each parameter    
    Returns
    ----------
    parameter: dictionary
    """

    if len(name) != len(minimum) or len(name) != len(maximum) or len(name) != len(value):
        print('Some limits are missing. %d parameteers with %d values, %d minima and %d maxima'
            %(len(name),len(value),len(minimum),len(maximum)))
        return None

    parameters = {}
    for par_name,par_value,par_min,par_max in zip(name,value,minimum,maximum):
        parameters.update({par_name : {'value':par_value,'min':par_min,'max': par_max}})

    return parameters


def update_parameters_values(parameters,value):
    '''Used in the fitting routine to update the parameter values in the correct order'''
    for par_name,par_value in zip(parameters.keys(),value):
        parameters[par_name]['value'] = par_value

    return parameters


def bin_data(data,binning_map):
    '''Used in the fitting routine to bin the data. Returns the data if no binning is present'''

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
    '''Used in the fitting routine to mask the data. Returns the data if no mask is present'''

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