import numpy as np
import os.path

import matplotlib.pylab as plt

from astropy.io import fits 
from astropy.cosmology import WMAP9 as cosmo
from reproject import reproject_interp
from astropy.convolution import convolve,Gaussian2DKernel

__all__ = ['Observation']

class Observation(object):
    """Data to be fitted.

    This class handles the data that is going to be fitted with the models. It includes
    the image (data) to be fit, the associated uncertainty (optionally ), the redshift 
    of the object and the seeing.

    Parameters
    ----------
    data_path: str
        Path to the fits file containing the data.
    unc_path: str
        Path to the fits file containing the associated uncertainty (optional).
    z: float 
        Redshift.
    seeing: float
        Seeing (in pixels). To be used in the model convolution.
    data: float array
        Data to be fitted. It is read from the data_path
    unc: float array
        Data uncertainted. It is read from the unc_path. 

    """
 

    def __init__(self,data_path,z,unc_path=None,seeing=0):

        self.z =  z
        self.seeing = seeing

        if os.path.isfile(data_path):
            self.data_path = data_path
        else:
            print('Data file %s not found'%data_path)
        if unc_path is not None:
            if os.path.isfile(data_path):
                self.unc_path  = unc_path
            else:
                print('Uncertainty file %s not found'%data_path)
        else:
            self.unc_path  = unc_path  

        # Load data and header
        self.data = fits.getdata(data_path)
        self.header = fits.getheader(data_path)

        # Load uncertainty if available. Assume all ones if not
        if unc_path is not None:
            self.unc = fits.getdata(unc_path)
        else:
            self.unc = np.ones_like(self.data)


    def info(self):
        """ Prints the data file, redshift and seeing."""

        print('Data: %s' %self.data_path)
        print('Redshift: %s' %self.z)
        print('Seeing (in pixels): %s' %self.seeing)

        if self.unc_path is None:
            print('No uncertainty image found. Using uniform uncertainty.')
        else:
            print('Uncertainty: %s'%self.unc_path)

        # Complete with binning and mask

    def plot(self,data_lim=None,unc_lim=None):
        """ Plots the data and uncertainty.
        Parameters
        ----------
        data_lim: (float,float)
            minimum and maximum values of the data colour bar. 
        unc_lim: (float,float)
            minimum and maximum values of the uncertainty colour bar. 
        """

        fig, ax = plt.subplots(1,2,figsize=(10,3))
        ax[0].set_title('Data')
        ax[1].set_title('Uncertainty')
        if data_lim is None:
            cax0 = ax[0].imshow(self.data,origin='lower')
        else:
            cax0 = ax[0].imshow(self.data,origin='lower',vmin=data_lim[0],vmax=data_lim[1])
        if unc_lim is None:
            cax1 = ax[1].imshow(self.unc,origin='lower')
        else:
            cax1 = ax[1].imshow(self.unc,origin='lower',vmin=unc_lim[0],vmax=unc_lim[1])
        plt.colorbar(cax0,ax=ax[0],fraction=0.03)
        plt.colorbar(cax1,ax=ax[1],fraction=0.03)

