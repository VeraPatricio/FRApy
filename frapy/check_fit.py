import numpy as np
import matplotlib.pylab as plt
import pickle
import corner

from .utils import bin_data,mask_data



__all__ = ['Output']

class Output(object):
    """ Allows the output of *fit_model* to be inspected
    
    Reads the pickle output and allows to plot:

        . the walkers positions at each iteration to check for convergence
        . a corner plot of the results
        . the 50th, 16th and 84th percentiles (mean and +/- 1 sigma)
 
    Parameters
    ----------
    outfile: str
        The name of the pickle file being inspected (without the '.pickle' extension)
    """

    def __init__(self,outfile):

        self.outfile = outfile
    
        results = pickle.load(open(self.outfile+".pickle",'rb'))    

        self.chain = results['chain']
        self.lnprobability = results['lnprobability'] 
        self.parameters_names = results['parameters_order']
        self.obs = results['obs']
        self.model = results['model']
        self.binning_map = results['bin_map']
        self.mask = results['mask']
        self.ndim = len(results['parameters_order'])

    def check_convergence(self):
        ''' Plots the walkers positions at each iteration for each parameter as well as 
        the value of the log-likelihood probability for each iteration.
        '''

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
        """ Makes a corner plot of the results. 
        Only uses iterations after 'start' 
        """       
        samples = self.chain[:, start:, :].reshape((-1, self.ndim))

        #best_par = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))

        fig = corner.corner(samples, labels=self.parameters_names,quantiles=[0.16, 0.5, 0.84], show_titles=True)


    def best_parameters(self,start=0):
        """Calculates the  16th, 50th and 84th percentiles for each parameter. 
        Only uses iterations after 'start' """

        samples = self.chain[:, start:, :].reshape((-1, self.ndim))

        best_par = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))

        best_parameters = {}
        for par_name,v in zip(self.parameters_names,best_par):
            best_parameters.update({par_name : {'value':v[0],'min':v[1],'max': v[2]}})
            print('%s %0.4f$^{+%0.4f}_{-%0.4f}$'%(par_name,v[0],v[1],v[2]))

        return best_parameters


    def plot_solution(self,best_parameters):
        """Given a dictionary with parameter names and values, plots the model.
 
        Parameters
        ----------
        best_parameters: dictionary
            Dictionary in the shape {parameter_name1:{value:X,min:Z,max:Z},parameter_name2:{value:X,min:Z,max:Z}}
            (from the check_fit.best_parameters function, for example).

        Returns
        ----------
        model: array float
            the model with the best parameters
        residuals: array float
            the residuals (data - model)
        """
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
        """Given a dictionary with parameter names and values, calculates the
        chi2, reduced chi2 (chi2/dof), the log-likelihood probability and the
        Bayesian Information Criteria (BIC) for the model with those parameters
        values.
 
        Parameters
        ----------
        best_parameters: dictionary
            Dictionary in the shape {parameter_name1:{value:X,min:Z,max:Z},parameter_name2:{value:X,min:Z,max:Z}}
            (from the check_fit.best_parameters function, for example).

        Returns
        ----------
        chi2/dof : float
            Reduced chi2
        """
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

        chi2_image = (dbata_obs-convolved_model)**2/data_unc**2
        chi2 = np.nansum(chi2_image[np.where(np.isfinite(chi2_image))])
        inv_sigma2 = 1.0/(data_unc**2)
        lnp = -0.5*(np.nansum((data_obs-convolved_model)**2*inv_sigma2- np.log(inv_sigma2*np.sqrt(2*np.pi))))
        bic = free_par*np.log(len(data_obs)) - 2*lnp

        print('Chi2: %0.2f'%chi2)
        print('Chi2/dof: %0.2f'%(chi2/(len(data_obs)-free_par)))
        print('Loglikelihood: %d'%lnp)
        print('BIC: %d'%bic)
        
        return chi2/(len(data_obs)-free_par)
