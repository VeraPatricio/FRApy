"""
The BaseModel handles all the lensing part, producing a distance map 
that is used by all the other models. At the moment, the following
models are implemented:

    * Metallcity gradient

    * TO BE COMPLETED
"""

import numpy as np
import os.path

import matplotlib.pylab as plt

from astropy.io import fits 
from astropy.cosmology import WMAP9 as cosmo
from reproject import reproject_interp
from astropy.convolution import convolve,Gaussian2DKernel

__all__ = ['Metallicity_Gradient']

class BaseModel(object):
    """ Global lensing model to be used in all other Model classes .

    This class prepares the displacement maps to be used with a particular object (i.e. 
    at a particular redshift) and observations (i.e. aligns the maps with the data).

    The main output is a distance map, in kiloparsecs, and an azimuthal map that serve
    as base for all the models being fit (metallicity gradient, velocity...)

    Parameters
    ----------
    z_lens: float
        Redshift of the gravitational lens.
    dplx_path: str
        Path to the fits file containing displacements in the x direction.
    dply_path: str
        Path to the fits file containing displacements in the y direction
    cx: int
        x position of the centre (in pixels)
    cy: int
        y position of the centre (in pixels)
    q: float
        axis ratio (a/b)
    pa: float
        Position angle (0 North, +90 East )
    lens_x: float array
        Lensing model (displacement in x direction) to be used to a particular object. 
        Created with the 'create_displacement_maps_for_object' method.
    lens_y: float array
        Lensing model (displacement in y direction) to be used to a particular object. 
        Created with the 'create_displacement_maps_for_object' method. 
    data: array
        An array with a realisation of a model made from the current parameter values.
    """

    def __init__(self,zlens,dplx_path,dply_path,cx=0,cy=0,q=1,pa=0):

        self.zlens = zlens
        if os.path.isfile(dplx_path):
            self.dplx_path = dplx_path
        else: 
            print('Displacement map %s not found'%dplx_path)
        if os.path.isfile(dply_path):
            self.dply_path = dply_path
        else: 
            print('Displacement map %s not found'%dply_path)
        self.cx = cx
        self.cy = cy
        self.q  = q
        self.pa = pa
        self.lens_x = None
        self.lens_y = None
        self.data = None


    def lensing_info(self):
        """ Prints the lens redshift and displacement maps origin"""    
        print('Lens redshift: %0.4f'%self.zlens)
        print('Displacement map (x): %s' %self.dplx_path)
        print('Displacement map (y): %s' %self.dply_path)
            
               
    def create_displacement_maps_for_object(self,Observation,correct_z=True):
        """ Takes the more global displacement maps produced by a graviatational 
        lensing fitting code, and converts these maps to displacements at the redhisft 
        of the source being analysed and aligns and regrids the maps to the data.
        The lens_x and lens_y attributes are created with this function.
        """
        # Load displacement maps
        dplx = fits.getdata(self.dplx_path)
        dply = fits.getdata(self.dply_path)
        dpl_header = fits.getheader(self.dplx_path)

        # Normalise it to the redshift of the object
        if correct_z:
            dls = cosmo.angular_diameter_distance_z1z2(self.zlens,Observation.z)
            ds = cosmo.angular_diameter_distance(Observation.z)
            dplx *= (dls/ds).value
            dply *= (dls/ds).value

        # Convert displacement maps to arcseconds 
        step = abs(dpl_header['CDELT2'])*3600.0
        x,y = np.meshgrid(xrange(dplx.shape[0]),xrange(dplx.shape[1]))
        sx =  x*step - dplx
        sy =  y*step - dply  

        # Re-grid the displacement map to the data format and to kpc
        data_header = fits.getheader(Observation.data_path)
        lens_x, _ = reproject_interp((sx,dpl_header),data_header) 
        lens_y, _ = reproject_interp((sy,dpl_header),data_header)
        lens_x = lens_x/cosmo.arcsec_per_kpc_proper(Observation.z).value
        lens_y = lens_y/cosmo.arcsec_per_kpc_proper(Observation.z).value

        self.lens_x = lens_x
        self.lens_y = lens_y


    def convolve_with_seeing(self,seeing):
        """Convolves a model with a Gaussian with width (sigma) 'seeing'."""
        return convolve(self.data, Gaussian2DKernel(stddev=(seeing)),boundary='extend')


    def make_distance_map(self):
        """Produces a distance map, in kpc, centrered in 'cx','cy' and assuming a ratio of 
        'q' between the minor and major axis, with the major axis in the 'pa' direction. """

        if np.all(self.lens_x) == None or np.all(self.lens_y) == None:
            print('No displacement maps for a particular redshift were found.') 
            print('It is not possible to create a distance map without them.')
            print('Use the "create_displacement_maps_for_object" method first.')
            pass

        else:

            # Centre distance map
            cx = int(np.round(self.cx,0))
            cy = int(np.round(self.cy,0))
            center_x = self.lens_x[cy,cx] 
            center_y = self.lens_y[cy,cx]  
            centered_lens_x = self.lens_x - center_x
            centered_lens_y = self.lens_y - center_y

            # Create 2D model
            y_ip,x_ip = np.mgrid[:centered_lens_x.shape[0], : centered_lens_x.shape[1]]
            x_out = centered_lens_x[y_ip,x_ip]
            y_out = centered_lens_y[y_ip,x_ip]
            
            # Account for rotation (counter-clockwise, North 0, East 90)
            pa_rad = np.deg2rad(self.pa - 90.)
            x_rot = x_out*np.cos(pa_rad)+y_out*np.sin(pa_rad) 
            y_rot = y_out*np.cos(pa_rad)-x_out*np.sin(pa_rad)
            dist = np.sqrt((x_rot)**2+((y_rot)/self.q)**2)
                    
            return dist

    def make_azimuthal_map(self):
        print('To be done')


    def plot(self):
        """Plots the model"""
        if np.all(self.data) == None:
            print('No model has been created yet. Use "make_model" method before "plot".')
        else:
            fig, ax = plt.subplots(1,1,figsize=(5,3))
            ax.set_title('Model')
            cax0 = ax.imshow(self.data,origin='lower')
            plt.colorbar(cax0,ax=ax,fraction=0.03)


class Metallicity_Gradient(BaseModel):
    """ Linear metallicity gradient.  

    This model inherits the distance maps attributes (cx,cy,q and pa), from which the metallicity
    at each point is calculated assuming a gradient and a central metallicity value:

    Z(r) = Delta Z * r + Z_0 

    with r the radius in kpc, Delta Z the gradient in dex/kpc, Z_0 the central metallicity. 
    
    Parameters
    ----------
    cx: int
        x position of the centre (in pixels)
    cy: int
        y position of the centre (in pixels)
    q: float
        axis ratio (a/b)
    pa: float
        Position angle (0 North, +90 East )
    z_grad: float
        Gradient in dex/kpc.
    z_0: float
        Central metallicity value (value at cx,cy)
    """

    def __init__(self,zlens,dplx_path,dply_path,cx=0,cy=0,q=1,pa=0,z_grad = -1,z_0 = 0):
        BaseModel.__init__(self,zlens,dplx_path,dply_path,cx=0,cy=0,q=1,pa=0)
        self.z_grad = z_grad
        self.z_0 = z_0

    def model_name():
        """Returns the model's name"""
        return 'metallicity_gradient'

    def model_parameters(self,verbose=True):
        """Returns the model's parameters"""
        if verbose:
            print('cx: x position of the centre (in pixels)')
            print('cy: y position of the centre (in pixels)')
            print('q: axis ratio (a/b)')
            print('pa: position angle (in degrees)')
            print('z_grad : gradient in dex/kpc')
            print('z_0: central metallicity')
        return ['cx','cy','q','pa','z_grad','z_0']

    def print_parameter_values(self):
        """Returns the model's parameters values"""
        print('cx: %d'%self.cx)
        print('cy: %d'%self.cy)
        print('q: %0.2f'%self.q)
        print('pa: %0.2f'%self.pa)
        print('z_grad: %0.2f'%self.z_grad)
        print('z_0: %0.2f'%self.z_0)


    def make_model(self):
        """ Makes a model using the current parameters' values and stores it 
        in the 'data' attribute"""
        distance_map = self.make_distance_map()
        grad = distance_map * self.z_grad + self.z_0
        self.data = distance_map * self.z_grad + self.z_0
        return grad


    def update_model_parameters(self,par):
        """Updates the parameters of the model.

        Parameters
        ----------
        par: dictionary
            dictionary in the shape {'name':parameter_name, 'value':parameter value}
        """
        for name in par.keys():
            if name == 'cx':
                self.cx = par[name]['value']
            if name == 'cy':
                self.cy = par[name]['value']
            if name == 'q':
                self.q= par[name]['value']
            if name == 'pa':
                self.pa = par[name]['value']
            if name == 'z_grad':
                self.z_grad = par[name]['value']
            if name == 'z_0':
                self.z_0 = par[name]['value']

class Metallicity_Gradient_with_Flatenning(BaseModel):
    """ Linear metallicity gradient with a break at the outter radius.  

    This model inherits the distance maps attributes (cx,cy,q and pa), from which the metallicity
    at each point is calculated assuming a gradient, a central metallicity value and a radius at
    which the flatenning happens:


    Z(r) = Delta Z * r + Z_0, for r <r_flat
         = cte, for r <r_flat

    with r the radius in kpc, Delta Z the gradient in dex/kpc, Z_0 the central metallicity and
    z_flat the radius at which the flatenning occurs.  
    
    Parameters
    ----------
    cx: int
        x position of the centre (in pixels)
    cy: int
        y position of the centre (in pixels)
    q: float
        axis ratio (a/b)
    pa: float
        Position angle (0 North, +90 East )
    z_grad: float
        Gradient in dex/kpc.
    z_0: float
        Central metallicity value (value at cx,cy)
    """

    def __init__(self,zlens,dplx_path,dply_path,cx=0,cy=0,q=1,pa=0,z_grad = 1,z_0 = 0):

        Model.__init__(self,zlens,dplx_path,dply_path,cx=0,cy=0,q=1,pa=0)
        self.z_grad = z_grad
        self.z_0 = z_0
        self.r_flat = r_flat

    def model_name():
        return 'metallicity_gradient'

    def make_model(self):
        distance_map = make_distance_map(self)
        grad = distance_map * self.z_grad + self.z_0
        outter_part = np.where(distance_map > self.r_flat)
        anulli = np.where((distance_map > 0.8*self.r_flat) & (distance_map < 1.2*self.r_flat))
        grad[outter_part] = np.mean(grad[anulli])
        return grad
