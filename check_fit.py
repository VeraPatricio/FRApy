"""
This module contains the Output class, to check the output of the fit.
"""

import numpy as np
import matplotlib.pylab as plt
import pickle
import corner

### Re-write this
from fit_model import bin_data,mask_data
from obs_and_model import Observation,Model,Metallicity_Gradient

class Output(object):

    def __init__(self,outfile):

        self.outfile = outfile
    
        results = pickle.load(open(self.outfile+".pickle",'rb'))    

        self.chain = results['chain']
        self.lnprobability = results['lnprobability'] 
        self.parameters = results['parameters']
        self.parameters_names = results['parameters_order']
        self.obs = results['obs']
        self.model = results['model']
        self.binning_map = results['bin_map']
        self.mask = results['mask']
        self.ndim = len(results['parameters'].keys())

    def check_convergence(self):

        nwalk = self.lnprobability.shape[0]

        fig1, ax = plt.subplots(1,self.ndim+1,figsize=(14,4))
        ax = ax.ravel()
        for j in range(nwalk):
            ax[0].plot(self.lnprobability[j, :])
            ax[0].set_title('lnP')
        
        for i in range(self.ndim):
            for j in range(nwalk):
                ax[i+1].plot(self.chain[j, :, i])
                ax[i+1].set_title(self.parameters_names[i])


    def make_cornerplot(self,start=0):
        """ start: chain steps below this value are not included in the plot"""
       
        samples = self.chain[:, start:, :].reshape((-1, self.ndim))

        best_par = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))

        fig = corner.corner(samples, labels=self.parameters_names,truths=np.array(best_par).T[0],\
                            quantiles=[0.16, 0.5, 0.84], show_titles=True)


    def best_parameters(self,start=0):

        samples = self.chain[:, start:, :].reshape((-1, self.ndim))

        best_par = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))

        best_parameters = self.parameters.copy()
        for par_name,v in zip(self.parameters_names,best_par):
            best_parameters.update({par_name : {'value':v[0],'min':v[1],'max': v[2]}})
            print('%s %0.4f$^{+%0.4f}_{-%0.4f}$'%(par_name,v[0],v[1],v[2]))

        return best_parameters


    def plot_solution(self,best_parameters):

        # Mask data if mask is available
        data_obs, data_unc = mask_data(self.obs,self.mask)

        # Update model
        model = self.model
        model.update_model_parameters(best_parameters)
        dummy_model = model.make_model()
        convolved_model = model.convolve_with_seeing(seeing=self.obs.seeing/2.355)
        residuals = self.obs.data - convolved_model

        fig, ax = plt.subplots(1,4,figsize=(10,5))
        ax[0].set_title('Data')
        ax[1].set_title('Model')
        ax[2].set_title('Convolved Model')
        ax[3].set_title('Residuals\n(Data-ConvolvedModel)')

        ax[0].imshow(self.obs.data,origin='lower')
        ax[1].imshow(model.data,origin='lower')
        ax[2].imshow(convolved_model,origin='lower')
        cax = ax[3].imshow(residuals,origin='lower',cmap='PiYG')
        plt.colorbar(cax,ax=ax[3],fraction=0.03)            

        return model, residuals


    def goodness_of_fit(self,best_parameters):

        free_par = len(self.parameters_names)

        # Mask data if mask is available
        data_obs, data_unc = mask_data(self.obs,self.mask)

        # Update model
        model = self.model
        model.update_model_parameters(best_parameters)
        dummy_model = model.make_model()
        convolved_model = model.convolve_with_seeing(seeing=self.obs.seeing/2.355)
        residuals = self.obs.data - convolved_model

        # If there is binning, bin data and model 
        data_obs = bin_data(self.obs.data,self.binning_map)
        data_unc = bin_data(self.obs.unc,self.binning_map)
        convolved_model = bin_data(convolved_model,self.binning_map)

        chi2 = np.nansum((data_obs-convolved_model)**2/data_unc**2)
        inv_sigma2 = 1.0/(data_unc**2)
        lnp = -0.5*(np.nansum((data_obs-convolved_model)**2*inv_sigma2- np.log(inv_sigma2*np.sqrt(2*np.pi))))
        bic = free_par*np.log(len(data_obs)) - 2*lnp

        print('Chi2: %0.2f'%chi2)
        print('Chi2/dof: %0.2f'%(chi2/(len(data_obs)-free_par)))
        print('Loglikelihood: %d'%lnp)
        print('BIC: %d'%bic)
        
        return chi2/(len(data_obs)-free_par)