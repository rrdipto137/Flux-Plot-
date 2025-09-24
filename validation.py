import openmc
import numpy as np

# --- materials ---
# Example: Low-enriched uranium metal (LEU) with 4.95% U-235, trace U-234
# Adjust fractions to match your approved, non-weaponizable material.
u_leu = openmc.Material(name='LEU metal (example)')
#u_leu.add_nuclide('U234', 0.0102)   # trace, adjust or set to 0 if unknown
u_leu.add_nuclide('U235', 90)   # 4.95% enrichment (example)
u_leu.add_nuclide('U238', 10)     # remainder
u_leu.set_density('g/cm3', 18.7398)    # uranium metal density ~19 g/cc

materials = openmc.Materials([u_leu])
materials.export_to_xml()

# --- geometry ---
sphere_radius =7# cm (requested geometry size, safe material)
world_size = 100.0    # cm half-length

# Surfaces
sphere = openmc.Sphere(r=sphere_radius)

left   = openmc.XPlane(x0=-world_size, boundary_type='vacuum')
right  = openmc.XPlane(x0= world_size, boundary_type='vacuum')
bottom = openmc.YPlane(y0=-world_size, boundary_type='vacuum')
top    = openmc.YPlane(y0= world_size, boundary_type='vacuum')
front  = openmc.ZPlane(z0=-world_size, boundary_type='vacuum')
back   = openmc.ZPlane(z0= world_size, boundary_type='vacuum')

# Regions
inside_sphere  = -sphere
outside_sphere = +sphere
world_box = +left & -right & +bottom & -top & +front & -back

# Cells
fuel_cell  = openmc.Cell(name='LEU_Sphere', fill=u_leu, region=inside_sphere)
void_cell  = openmc.Cell(name='Void', region=outside_sphere & world_box)

# Root universe / geometry
root_universe = openmc.Universe(cells=[fuel_cell, void_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# --- eigenvalue run settings ---


src = openmc.Source()
src.space = openmc.stats.Point((0.0, 0.0, 0.0))
src.angle = openmc.stats.Isotropic()

src.particle = 'neutron'

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.source = src
settings.batches = 300
settings.inactive = 50
settings.particles = 1000000
settings.export_to_xml()

# --- run ---
openmc.run()

