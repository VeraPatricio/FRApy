"""
The BaseModel handles all the lensing part, producing a distance map 
that is used by all the other models. 
"""

import numpy as np
import os.path

import matplotlib.pylab as plt

from astropy.io import fits 
from astropy.cosmology import WMAP9 as cosmo
from reproject import reproject_interp
from astropy.convolution import convolve,Gaussian2DKernel

__all__ = ['BaseModel','Metallicity_Gradient','Metallicity_Gradient_Constant_Centre','Velocity_Arctangent']

class BaseModel(object):
    """ Global lensing model to be used in all other Model classes .

    This class prepares the deflection maps to be used with a particular object (i.e. 
    at a particular redshift) and observations (i.e. aligns the maps with the data).

    The main output is a distance map, in kiloparsecs, and an azimuthal map that serve
    as base for all the models being fit (metallicity gradient, velocity...)

    Parameters
    ----------
    z_lens: float
        Redshift of the gravitational lens.
    dfx_path: str
        Path to the fits file with the x deflection.
    dfy_path: str
        Path to the fits file with the y deflection.
    cx: int
        x position of the centre (in pixels)
    cy: int
        y position of the centre (in pixels)
    q: float
        axis ratio (b/a)
    pa: float
        Position angle (0 North, +90 East )
    df_ang: float
        Angle between x axis and North (measured anti-clockwise) in the deflection maps.
    project_x: float array
        Lensing model (deflection in x direction) to be used to a particular object. 
        Created with the 'create_deflection_maps_for_object' method.
    project_y: float array
        Lensing model (deflection in y direction) to be used to a particular object. 
        Created with the 'create_deflection_maps_for_object' method. 
    data: array
        An array with a realisation of a model made from the current parameter values.
    conv_data: array
        An array with a realisation of a model made from the current parameter values,
        convolved by the seeing of observations.
    """

    def __init__(self,zlens,dfx_path,dfy_path,df_ang=0,cx=0,cy=0,q=1,pa=0):

        self.zlens = zlens
        if os.path.isfile(dfx_path):
            self.dfx_path = dfx_path
        else: 
            print('Deflection map %s not found'%dfx_path)
        if os.path.isfile(dfy_path):
            self.dfy_path = dfy_path
        else: 
            print('Deflection map %s not found'%dfy_path)
        self.cx = cx
        self.cy = cy
        self.q  = q
        self.pa = pa
        self.project_x = None
        self.project_y = None
        self.data = None
        self.conv_data = None
        self.df_ang = df_ang

    def lensing_info(self):
        """ Prints the lens redshift and deflection maps origin"""    
        print('Lens redshift: %0.4f'%self.zlens)
        print('Deflection map (x): %s' %self.dfx_path)
        print('dDflection map (y): %s' %self.dfy_path)
            
               
    def create_projection_maps(self,Observation,correct_z=True):
        """ Takes the more global deflection maps produced by a graviatational 
        lensing fitting code, and converts these maps to 'projection' maps, that
        maps where a pixel in source plane should be 'projected' in image plane,
        for this particular Observation. The project_x and project_y attributes 
        are created with this function.
        """

        # Load deflection maps
        dplx = fits.getdata(self.dfx_path)
        dply = fits.getdata(self.dfy_path)
        dpl_header = fits.getheader(self.dfx_path)

        # Normalise it to the redshift of the object
        if correct_z:
            dls = cosmo.angular_diameter_distance_z1z2(self.zlens,Observation.z)
            ds = cosmo.angular_diameter_distance(Observation.z)
            dplx *= (dls/ds).value
            dply *= (dls/ds).value

        # Convert deflection maps to arcseconds 
        step = abs(dpl_header['CDELT2'])*3600.0
        x,y = np.meshgrid(np.arange(dplx.shape[0]),np.arange(dplx.shape[1]))
        sx =  x*step - dplx
        sy =  y*step - dply  

        # Re-grid the deflection map to the data format and to kpc
        data_header = fits.getheader(Observation.data_path)
        project_x, _ = reproject_interp((sx,dpl_header),data_header) 
        project_y, _ = reproject_interp((sy,dpl_header),data_header)
        project_x = project_x/cosmo.arcsec_per_kpc_proper(Observation.z).value
        project_y = project_y/cosmo.arcsec_per_kpc_proper(Observation.z).value

        self.project_x = project_x
        self.project_y = project_y


    def convolve_with_seeing(self,seeing):
        """Convolves a model with a Gaussian with width (sigma) 'seeing'."""
        conv_model = convolve(self.data, Gaussian2DKernel(stddev=(seeing)),boundary='extend')
        self.conv_data = conv_model
        return conv_model


    def make_distance_map(self):
        """Produces a distance map, in kpc, centrered in 'cx','cy' and assuming a ratio of 
        'q' between the minor and major axis, with the major axis in the 'pa' direction. """

        if np.all(self.project_x) == None or np.all(self.project_y) == None:
            print('No projection maps for a particular redshift were found.') 
            print('It is not possible to create a distance map without them.')
            print('Use the "create_projection_maps" method first.')
            pass

        else:

            # Centre distance map
            cx = int(np.round(self.cx,0))
            cy = int(np.round(self.cy,0))
            center_x = self.project_x[cy,cx] 
            center_y = self.project_y[cy,cx]  
            centered_project_x = self.project_x - center_x
            centered_project_y = self.project_y - center_y

            # Create 2D model
            y_ip,x_ip = np.mgrid[:centered_project_x.shape[0], : centered_project_x.shape[1]]
            x_out = centered_project_x[y_ip,x_ip]
            y_out = centered_project_y[y_ip,x_ip]
            
            # Account for rotation (counter-clockwise, North 0, East 90)
            pa_rad = np.deg2rad(self.df_ang + self.pa)
            x_rot = x_out*np.cos(pa_rad)+y_out*np.sin(pa_rad) 
            y_rot = y_out*np.cos(pa_rad)-x_out*np.sin(pa_rad)
            dist = np.sqrt((x_rot*self.q)**2+(y_rot)**2)

            return dist

    def make_azimuthal_map(self):        
        """Produces an azimuthal map, in kpc, centrered in 'cx','cy' and assuming a ratio of 
        'q' between the minor and major axis, with the major axis in the 'pa' direction. """

        if np.all(self.project_x) == None or np.all(self.project_y) == None:
            print('No deflection maps for a particular redshift were found.') 
            print('It is not possible to create a distance map without them.')
            print('Use the "create_deflection_maps_for_object" method first.')
            pass

        else:

            # Centre distance map
            cx = int(np.round(self.cx,0))
            cy = int(np.round(self.cy,0))
            center_x = self.project_x[cy,cx] 
            center_y = self.project_y[cy,cx]  
            centered_project_x = self.project_x - center_x
            centered_project_y = self.project_y - center_y

            # Create 2D model
            y_ip,x_ip = np.mgrid[:centered_project_x.shape[0], : centered_project_x.shape[1]]
            x_out = centered_project_x[y_ip,x_ip]
            y_out = centered_project_y[y_ip,x_ip]
            
            # Account for rotation (counter-clockwise, North 0, East 90)
            pa_rad = np.deg2rad(self.df_ang + self.pa)
            x_rot = x_out*np.cos(pa_rad)+y_out*np.sin(pa_rad) 
            y_rot = y_out*np.cos(pa_rad)-x_out*np.sin(pa_rad)
            dist = np.sqrt((x_rot*self.q)**2+(y_rot)**2)

            ratio_x = x_rot/dist
            ratio_y = y_rot/dist
            ang = np.rad2deg(np.arccos(ratio_x))
            neg_ratio = np.where(ratio_y > 0)
            ang[neg_ratio] = 360 -  np.rad2deg(np.arccos(ratio_x[neg_ratio]))

            return ang

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

    def __init__(self,zlens,dfx_path,dfy_path,df_ang=0,cx=0,cy=0,q=1,pa=0,z_grad = -1,z_0 = 0):
        BaseModel.__init__(self,zlens,dfx_path,dfy_path,cx=cx,cy=cy,q=q,pa=pa,df_ang=df_ang)
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
        print('df_ang: %0.2f'%self.df_ang)
        print('z_grad: %0.4f'%self.z_grad)
        print('z_0: %0.4f'%self.z_0)

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

    def make_model(self):
        """ Makes a model using the current parameters' values and stores it 
        in the 'data' attribute"""
        distance_map = self.make_distance_map()
        grad = distance_map * self.z_grad + self.z_0
        self.data = distance_map * self.z_grad + self.z_0
        return grad


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

    def __init__(self,zlens,dfx_path,dfy_path,cx=0,cy=0,q=1,pa=0,z_grad = 1,z_0 = 0):

        Model.__init__(self,zlens,dfx_path,dfy_path,cx=0,cy=0,q=1,pa=0)
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


class Metallicity_Gradient_Constant_Centre(BaseModel):
    """ Linear metallicity gradient with a flatenning of the centre at r_flat.  

    This model inherits the distance maps attributes (cx,cy,q and pa), from which the metallicity
    at each point is calculated assuming a gradient and a central metallicity value:

    Z(r) = Delta Z * r + Z_0 

    with r the radius in kpc, Delta Z the gradient in dex/kpc, Z_0 the central metallicity. 

    For r < r_flat, the metallicity is constant.
    
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
    r_flat: float
        Radius that delimits the central zone where the metallicity is flat
    """

    def __init__(self,zlens,dfx_path,dfy_path,cx=0,cy=0,q=1,pa=0,z_grad = -1,z_0 = 0,r_flat=0.5,z_grad_inner = -1):
        BaseModel.__init__(self,zlens,dfx_path,dfy_path,cx=0,cy=0,q=1,pa=0)
        self.z_grad = z_grad
        self.z_0 = z_0
        self.r_flat = r_flat
        self.z_grad_inner = z_grad_inner


    def model_name():
        """Returns the model's name"""
        return 'Metallicity_Gradient_constant_centre'

    def model_parameters(self,verbose=True):
        """Returns the model's parameters"""
        if verbose:
            print('cx: x position of the centre (in pixels)')
            print('cy: y position of the centre (in pixels)')
            print('q: axis ratio (a/b)')
            print('pa: position angle (in degrees)')
            print('z_grad : gradient in dex/kpc')
            print('z_0: central metallicity')
            print('r_flat: maximum radius of the central (flat) region')
            print('z_grad_inner: metallicity gradient of the inner part')
        return ['cx','cy','q','pa','z_grad','z_0','r_flat','z_grad_inner']

    def print_parameter_values(self):
        """Returns the model's parameters values"""
        print('cx: %d'%self.cx)
        print('cy: %d'%self.cy)
        print('q: %0.2f'%self.q)
        print('pa: %0.2f'%self.pa)
        print('z_grad: %0.4f'%self.z_grad)
        print('z_0: %0.4f'%self.z_0)
        print('r_flat: %0.2f'%self.r_flat)
        print('z_grad_inner: %0.4f'%self.z_grad_inner)

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
            if name == 'r_flat':
                self.r_flat = par[name]['value']
            if name == 'z_grad_inner':
                self.z_grad_inner = par[name]['value']

    def make_model(self):
        """ Makes a model using the current parameters' values and stores it 
        in the 'data' attribute"""
        distance_map = self.make_distance_map()
        grad = distance_map * self.z_grad + self.z_0
        grad_inner = distance_map * self.z_grad_inner + self.z_0
        inner_part = np.where(distance_map < self.r_flat)
        grad[inner_part] = grad_inner[inner_part]
        self.data = grad
        return grad


class Velocity_Arctangent(BaseModel):
    """ Exponential velocity model.  

    This model inherits the distance and azimuthal maps, from which an arctangent model
    of the velocity at each point is calculated assuming the following formulae:

    V(r) = v_t \\frac{2}{\pi} arctan (\\frac{2r}{r_t})

    with r the radius in kpc, v_t the terminal velocity and r_t the transition radius. 

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
    v_t: float
        Terminal velocity in km/s.
    r_t: float
        transition radius in kpc.    
    """

    def __init__(self,zlens,dfx_path,dfy_path,cx=0,cy=0,q=1,pa=0,v_t =100,r_t=10):
        BaseModel.__init__(self,zlens,dfx_path,dfy_path,cx=0,cy=0,q=1,pa=0)
        self.v_t = v_t
        self.r_t = r_t

    def model_name():
        """Returns the model's name"""
        return 'arctangent_velocity'

    def model_parameters(self,verbose=True):
        """Returns the model's parameters"""
        if verbose:
            print('cx: x position of the centre (in pixels)')
            print('cy: y position of the centre (in pixels)')
            print('q: axis ratio (a/b)')
            print('pa: position angle (in degrees)')
            print('v_t : terminal velocity')
            print('r_t: transition velocity')
        return ['cx','cy','q','pa','v_t','r_t']

    def print_parameter_values(self):
        """Returns the model's parameters values"""
        print('cx: %d'%self.cx)
        print('cy: %d'%self.cy)
        print('q: %0.2f'%self.q)
        print('pa: %0.2f'%self.pa)
        print('v_t: %0.2f'%self.v_t)
        print('r_t: %0.2f'%self.r_t)


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
            if name == 'r_t':
                self.r_t = par[name]['value']
            if name == 'v_t':
                self.v_t = par[name]['value']

    def make_model(self):
        """ Makes a model using the current parameters' values and stores it 
        in the 'data' attribute"""

        distance_map = self.make_distance_map()
        theta_map = self.make_azimuthal_map()

        vel = self.v_t * (2/np.pi) * np.arctan(2*distance_map/self.r_t) * (1 - self.q**2) * np.cos(np.deg2rad(theta_map))
        self.data = vel

        return vel