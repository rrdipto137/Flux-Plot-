import openmc
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Remove previous XML and HDF5 files
os.system('rm *.xml *.h5')

# Set cross-sections path
openmc.config['cross_sections'] = '/media/rrdipto/Dipto/Research/Nuclear_Data_Libraries/endfb80/endfb-viii.0-hdf5/cross_sections.xml'


u_leu = openmc.Material(name='LEU metal')
u_leu.add_nuclide('U235', 0.9)  # 4.95% enrichment
u_leu.add_nuclide('U238', 0.1)  # Remainder
u_leu.set_density('g/cm3', 18.74)  # Uranium metal density
u_leu.id = 1  # Explicit material ID

materials = openmc.Materials([u_leu])
materials.export_to_xml()

# --- Geometry ---
sphere_radius = 3.0  # cm
world_size = 10.0  # cm half-length (reduced for efficiency)

# Surfaces
sphere = openmc.Sphere(r=sphere_radius)
left = openmc.XPlane(x0=-world_size, boundary_type='vacuum')
right = openmc.XPlane(x0=world_size, boundary_type='vacuum')
bottom = openmc.YPlane(y0=-world_size, boundary_type='vacuum')
top = openmc.YPlane(y0=world_size, boundary_type='vacuum')
front = openmc.ZPlane(z0=-world_size, boundary_type='vacuum')
back = openmc.ZPlane(z0=world_size, boundary_type='vacuum')

# Regions
inside_sphere = -sphere
outside_sphere = +sphere & +left & -right & +bottom & -top & +front & -back

# Cells
fuel_cell = openmc.Cell(cell_id=1, name='LEU_Sphere', fill=u_leu, region=inside_sphere)
void_cell = openmc.Cell(cell_id=2, name='Void', region=outside_sphere)

# Root universe / geometry
root_universe = openmc.Universe(cells=[fuel_cell, void_cell])
root_universe.plot(width=(8,8),color_by='material',pixels=[400,400])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# --- Eigenvalue Run Settings ---
src = openmc.Source()
src.space = openmc.stats.Point((0.0, 0.0, 0.0))
src.angle = openmc.stats.Isotropic()
src.particle = 'neutron'

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.source = src
settings.batches = 100
settings.inactive = 5
settings.particles = 100000
settings.temperature = {'method': 'interpolation'}  # Enable temperature effects
settings.export_to_xml()



#energy flux
energies = np.logspace(np.log10(10e2), np.log10(10e8), 500)
e_filter = openmc.EnergyFilter(energies)
tally1 = openmc.Tally(name='energy flux')
tally1.id=1
tally1.filters=[e_filter]
tally1.scores=['flux']

tallies = openmc.Tallies([tally1])
tallies.export_to_xml()


#plot flux spectrum
sp = openmc.StatePoint('statepoint.100.h5')

t = sp.tallies[tally1.id]

#energy_centers = 0.5 * (energies[:-1] + energies[1:])
#plt.loglog(energy_centers, flux, , label='flux')

flux = t.mean.flatten() 
plt.loglog(energies[:-1], flux )
#plt.title('energy_flux_distribution')
plt.grid()
plt.xlabel('Energy[ev]')
plt.ylabel('Flux[n-cm/source]')
plt.legend()

# normalizing flux to get real flux value in n/cm^2.sec


volume = (4/3)*np.pi*3*3*3

assembly_power = 1           #set urself 
assembly_power_in_ev = assembly_power/1.6E-19
energy_per_fission = 200E6
neutron_souce = assembly_power_in_ev/energy_per_fission/2.43

flux = neutron_souce* t.mean.flatten() /volume
plt.loglog(energies[:-1], flux)
#plt.title('energy_flux_distribution')
plt.grid()
plt.xlabel('Energy[ev]')
plt.ylabel('Flux[n/cm^2.sec]')
plt.legend()


#spatial flux (delete the statepoint and Summary first)
model = openmc.Model(materials=materials,geometry=geometry, settings=settings)
d = 700 #
# Instantiate a tally Mesh
mesh = openmc.RegularMesh()
mesh.dimension = [d, d]
mesh.lower_left = [-4, -4]          #set urself
mesh.upper_right = [+4, +4]         #set urself

# Instantiate tally Filter
mesh_filter = openmc.MeshFilter(mesh)

# Creating Conversion Factor of Flux
# Creating a tally of fission-q-recoverable
H = openmc.Tally()
H.scores = ['heating-local']
H.filter = [mesh_filter]

# Instantiate the Tally
tally = openmc.Tally(name='mesh tally')
tally.filters = [mesh_filter]
tally.scores = ['flux']

tallies = openmc.Tallies([H, tally])

model.tallies = tallies

# Run OpenMC
statepoint_filename = model.run()


# Move the statepoint File
ce_spfile = './statepoint.100.h5'
#os.rename(statepoint_filename, ce_spfile)
# Move the Summary file
ce_sumfile = './summary.h5'
#os.rename('summary.h5', ce_sumfile)

# Load the statepoint file
sp = openmc.StatePoint(ce_spfile, autolink=False)

H_value = sp.tallies[H.id].get_pandas_dataframe()['mean'][0]

# ================================================================
# Calculating the Volume of Mesh
mesh_height = 60
mesh_area = (148 / d)**2          #set urself

volume = mesh_height * mesh_area

# Calculating the Factor to Normalize the Flux
power = 1
H_Jsrc = 1.602e-19 * H_value

f = power / H_Jsrc
flux_normalization_factor = f/volume

# ================================================================



# Load the summary file in its new location
su = openmc.Summary(ce_sumfile)
sp.link_with_summary(su)

# Get the OpenMC fission rate mesh tally data
ce_mesh_tally = sp.get_tally(name='mesh tally')
ce_flux = ce_mesh_tally.get_values(scores=['flux'])

# Reshape array to 2D for plotting
ce_flux.shape = mesh.dimension

# Normalize to the average pin power
ce_flux /= np.mean(ce_flux[ce_flux > 0.])

# Force zeros to be NaNs so their values are not included when matplotlib calculates
# the color scale
ce_flux[ce_flux == 0.] = np.nan


plt.figure(dpi=1200)
plt.imshow(ce_flux*flux_normalization_factor, interpolation='none', origin='lower', cmap = 'twilight_shifted')
plt.colorbar()
plt.title('Flux Distribution')
plt.show()



df = sp.tallies[tally.id].get_pandas_dataframe()

### Converting (490000,) norm_flux to 2d norm_flux with shape (700,700)

norm_flux = df['mean']*flux_normalization_factor
# Convert norm_flux to NumPy array if it's not already
norm_flux_array = norm_flux.to_numpy()  # From Pandas Series to NumPy array

# Reshape to 2D array
norm_flux_2d = norm_flux_array.reshape((d, d))

'''
'''
from matplotlib import colormaps

preferred_colormaps = [                        
 'twilight',
 'twilight_shifted',
 'turbo',
 'gist_stern',
 'magma_r',
 'inferno_r',
 'plasma_r',
 'viridis_r',
 'cividis_r',
 'twilight_r',
 'twilight_shifted_r',
 'turbo_r',
 'Wistia_r',
 'coolwarm_r',
 'copper_r',
 'cubehelix_r',
 'gist_stern_r',
 'gist_yarg_r',
 'seismic_r',
 ]                              #set urself

for colormap in preferred_colormaps:
    plt.figure(dpi=600)                    #set urself
    plt.imshow(norm_flux_2d, interpolation='none', origin='lower', cmap=colormap)
    plt.colorbar()
    plt.title('Flux Distribution')
    plt.show()
'''

plt.figure(dpi=600)
plt.imshow(norm_flux_2d, interpolation='none', origin='lower', cmap='twilight_shifted')
plt.colorbar(orientation='vertical', label='neutrons/$cm^2$-s')
plt.show()



plt.imshow(ce_flux*flux_normalization_factor, cmap='twilight_shifted', origin='lower')
plt.colorbar(label='Flux [n/cm²/s]')
plt.title("Flux Distribution")
'''


