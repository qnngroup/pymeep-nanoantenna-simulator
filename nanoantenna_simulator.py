#Headers
import meep as mp
import numpy as np
import h5py
import yaml
import sys
from mpi4py import MPI
from scipy.io import loadmat
from scipy.io import savemat
import meep.materials as meepmat

#Constants
# -- Units noted.  
# -- Fundamental unit length of 1-um assumed for meep calculations (see meep manual)
tcon = 3.33564 #time conversion to meep units (fs/meep time)
q =  1.6e-19 # electron charge in coulomb
m =  9.10938356e-31 # electron mass in kg 

#Functions for calculating index of refraction -- meep style
def calc_sig_d(n, k, fcen):
    eps = (n + 1j*k)**2
    eps_r = np.real(eps)
    return 2*np.pi*fcen*np.imag(eps)/eps_r

def calc_eps_r(n, k):
    eps = (n + 1j*k)**2
    return np.real(eps)

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def run_simulation(save_prefix, reference=False):

    #Open the settings yaml file
    with open(save_prefix + '_settings.yaml', 'r') as f:
        settings = yaml.safe_load(f)

    #----------------------------------
    # Load User Settings
    #----------------------------------

    # -- Simulation Type Toggles
    want_avgeps_on = settings['want_avgeps_on']
    want_visualize_on = settings['want_visualize_on']

    time_sample_rate = np.double(settings['time_sample_rate'])
    fcen = np.double(settings['fcen'])
    df = np.double(settings['df'])
    nfreq = int(settings['nfreq'])

    dpml = np.double(settings['dpml'])
    sx = np.double(settings['sx'])
    sy = np.double(settings['sy'])
    sz = np.double(settings['sz'])
    res = np.double(settings['res'])

    metal_thickness = np.double(settings['metal_thickness'])
    altitude = np.double(settings['altitude'])
    r_curvature = np.double(settings['r_curvature'])
    
    substrate_material_name = settings['substrate_material_name']
    metal_material_name = settings['metal_material_name']


    #---------------------------
    # FDTD Simulation Setup
    #---------------------------

    want_structure_on = not(reference)
    
    substrate_material = []
    metal_material = []

    air = mp.Medium(epsilon=1.0)
    
    if want_visualize_on:
        substrate_material=mp.Medium(epsilon=2.0)
        metal_material=mp.Medium(epsilon=3.0)

        save_prefix = save_prefix + "-visualize"
        
    else:
            
        #Define materials from meepmat database:
        substrate_material = eval('meepmat.' + substrate_material_name)
        metal_material = eval('meepmat.' + metal_material_name)

    # -- Define the Simulation Cell -- 
    #Calculate actual cell size after accounting for pmls
    # (only in Z here... assumed periodic in other direction)
    sZ = 2*dpml + sz
    sX = sx
    sY = sy

    cell_size = mp.Vector3(sX, sY, sZ)
       
    # Calculate convenience constants for setting the goemetry:
    base = altitude*0.75
    theta = np.arctan(8/3)
    r_top = r_curvature*np.tan(theta)
    height_cutout = r_top*np.cos(np.pi/2.0 - theta)
    half_base_cutout = r_top*np.sin(np.pi/2 - theta)
    altitude_adjusted = altitude - height_cutout
    side = np.sqrt((0.5*base - half_base_cutout)**2.0 + altitude_adjusted**2.0)
    cutout_center_offset = half_base_cutout + np.cos(theta)*side/2.0

    # -- Set Geometry --
    geometry = []

    if(want_structure_on):

        #Triangle
        geometry.append(mp.Block(center=mp.Vector3(0, 0, -1*metal_thickness/2.0),
                                 size=mp.Vector3(base,
                                                 altitude_adjusted,
                                                 metal_thickness),
                                 material=metal_material))

        geometry.append(mp.Block(center=mp.Vector3(base/4.0 + cutout_center_offset,
                                                   0,
                                                   -1*metal_thickness/2.0),
                                 e2=mp.Vector3(-1*base/2.0, altitude, 0),
                                 size=mp.Vector3(base/2.0, side, metal_thickness),
                                 material=air))

        geometry.append(mp.Block(center=mp.Vector3(-1*base/4.0 - cutout_center_offset,
                                               0,
                                               -1*metal_thickness/2.0),
                             e2=mp.Vector3(base/2.0, altitude, 0),
                             size=(base/2.0, side, metal_thickness),
                             material=air))

        #Round the triangle's tip
        geometry.append(mp.Cylinder(center=mp.Vector3(0,
                                                  altitude_adjusted/2.0 - r_curvature*np.cos(theta),
                                                  -1*metal_thickness/2.0),
                                radius=r_curvature,
                                height=metal_thickness,
                                material=metal_material));

        #Substrate
        geometry.append(mp.Block(center=mp.Vector3(0, 0, sZ*0.25),
                             size=mp.Vector3(mp.inf, mp.inf, sZ*0.5),
                             material=substrate_material))
        

    else:

        geometry.append(mp.Block(center=mp.Vector3(0, 0, 0),
                                 size=mp.Vector3(mp.inf, mp.inf, mp.inf),
                                 material=air))

    # -- Periodic Bdr. Cond -- 
    k_point = mp.Vector3(0, 0, 0)

    # -- Symmetries --
    symmetries = mp.Mirror(mp.X_DIR)

    # -- PML Layers --
    pml_layers = [mp.PML(thickness=np.double(dpml), direction=mp.Z)]


    # -- Near-Field Monitor --
    near_field_monitor = mp.Volume(center=mp.Vector3(0, 0, -1*metal_thickness/2.0),
                                   size=mp.Vector3(sx, sy, 0))


    # -- Define Sources --

    #Y-Polarized Gaussian source next to the dpml at -0.5 sz
    sources = [mp.Source(mp.GaussianSource(frequency=fcen, fwidth=df),
                         component=mp.Ey,
                         center=mp.Vector3(0, 0, -0.125*sz),
                         size=mp.Vector3(sx,sy,0))]

    # -- Load Settings Into Simulation -- 
    sim = mp.Simulation(resolution=res,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        k_point=k_point,
                        geometry=geometry,
                        dimensions=3,
                        filename_prefix=save_prefix,
                        eps_averaging=want_avgeps_on)

    # -- Transmitted Flux -- 
    trans_fr = mp.FluxRegion(center=mp.Vector3(0, 0, 0.25*sz),
                             size=mp.Vector3(sx, sy, 0),
                             direction=mp.Z)
    trans = sim.add_flux(fcen, df, nfreq, trans_fr)

    # -- Reflected Flux -- 
    refl_fr = mp.FluxRegion(center=mp.Vector3(0, 0, -0.25*sz),
                            size=mp.Vector3(sx, sy, 0),
                            direction=mp.Z)
    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    # -- Run Simulation --
    # - Different depending on whether this is a background,
    #   full simulation, or a visualization run
    
    if want_visualize_on:
        
        sim.run(mp.at_beginning(mp.output_epsilon),
                until=0)

    else:

        if want_structure_on:
            sim.load_minus_flux("refl-flux-ref", refl)

            sim.run(mp.at_beginning(mp.output_epsilon),
                    mp.in_volume(near_field_monitor,
                                 mp.to_appended("tEx_xy",
                                                mp.at_every(1.0/(fcen + 0.5*df)/time_sample_rate,
                                                            mp.output_efield_x))),
                    mp.in_volume(near_field_monitor,
                                 mp.to_appended("tEy_xy",
                                                mp.at_every(1.0/(fcen + 0.5*df)/time_sample_rate,
                                                            mp.output_efield_y))),
                    mp.in_volume(near_field_monitor,
                                 mp.to_appended("tEz_xy",
                                                mp.at_every(1.0/(fcen + 0.5*df)/time_sample_rate,
                                                            mp.output_efield_z))),

                    until=mp.stop_when_fields_decayed(dt=10, c=mp.Ey, pt=mp.Vector3(0, 0, sz*0.25), decay_by=1e-3))

            # -- Save output fluxes -- 
            flux_freqs = mp.get_flux_freqs(refl)
            refl_flux_out = mp.get_fluxes(refl)
            trans_flux_out = mp.get_fluxes(trans)

            savemat(save_prefix + '_output_flux.mat',
                    {'flux_freqs':flux_freqs,
                     'refl_flux_out':refl_flux_out,
                     'trans_flux_out':trans_flux_out})


        else:
            
            sim.run(mp.at_beginning(mp.output_epsilon),
                    mp.in_volume(near_field_monitor,
                                 mp.to_appended("bEy_xy",
                                                mp.at_every(1.0/(fcen + 0.5*df)/time_sample_rate,
                                                            mp.output_efield_y))),

                    until=mp.stop_when_fields_decayed(dt=10, c=mp.Ey, pt=mp.Vector3(0, 0, sz*0.25), decay_by=1e-3))

            sim.save_flux("refl-flux-ref", refl)                                 

# -- Command Line Argument Parser -- 
# - This part is essential for being able to run the code from the command line.
#   with appropriate settings flags.
# - It is needed for running the code with MPI capability to take advantage of
#   multiple thread operation.  
if __name__ == "__main__":

    if(len(sys.argv) == 3):

        if(sys.argv[2] == '-b'):

            run_simulation(sys.argv[1], True)

        else:

            print('Input ' + sys.argv[2] + ' not valid.')
            

    elif(len(sys.argv) == 2):

        run_simulation(sys.argv[1])

    else:

        print('Only maximum of two arguments allowed!')

        

        

