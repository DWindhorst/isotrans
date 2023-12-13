# loads libraries
import cmf
from datetime import datetime, timedelta
from math import exp
import numpy as np
from numpy.random import rand

from iso_time import iso_time
from iso_layers import iso_layer
from Boundary_conditions import BoundaryCondition, iso_atmosphere
import solve_iso_transport as iso
from Visualize import Visualize

import matplotlib.pyplot as plt

# Setup CMF
# A function to create a retention curve (used for the whole profile)
# for a specific depth, using an exponential decline function for Ksat with depth

# Likelihood Tests (Mathieu and Bariac (1996).)
ignorealphai = True
ignorealphaik = True
ignoredl = True
ignoredv = True

Braud_Ksat = 0.0106272  # in m/d
Braud_phi = 0.35  # porosity
Braud_alpha = 0.193  # Scale value of the water pressure(m)
Braud_n = 2.22
Braud_m = 0.099
Braud_theta_r = 0.01
Braud_eta = 9.14

VGM_retention_curve = cmf.VanGenuchtenMualem(Ksat=Braud_Ksat,
                                             phi=Braud_phi,
                                             alpha=Braud_alpha,
                                             n=Braud_n,
                                             m=Braud_m,
                                             theta_r=Braud_theta_r)

# Create a project
p = cmf.project()

# Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
c = p.NewCell(0, 0, 0, 1000)

# Customize cell
# Top layer thickness of e10-5 m as per SISPAT_iso
lower_boundaries_of_layer = [0.0001, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                             0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]


for d in lower_boundaries_of_layer:
    c.add_layer(d, VGM_retention_curve)
    c.saturated_depth = d

# Create a surfacewater storage
c.surfacewater_as_storage()
# use Richards connection
c.install_connection(cmf.Richards)

# Create a Neumann Boundary condition connected to top layer to account for flows due to evaporation
E_act = cmf.NeumannBoundary_create(c.layers[0])

# Create a List of Neumann Boundary condition and connect each one to a layer to account for possible lateral flow
lateral_flows = cmf.NeumannBoundary_list()
for i_layer in c.layers:
    i_lateral_flow = cmf.NeumannBoundary_create(i_layer)
    i_lateral_flow.flux = cmf.timeseries_from_scalar(0.0)
    lateral_flows.append(i_lateral_flow)

# Create a List of Neumann Boundary condition and connect each one to a layer to account for possible flow due to transpiration/root extraction
transpiration_flows = cmf.NeumannBoundary_list()
for i_layer in c.layers:
    i_transpiration_flow = cmf.NeumannBoundary_create(i_layer)
    i_transpiration_flow.flux = cmf.timeseries_from_scalar(0.0)
    lateral_flows.append(i_transpiration_flow)

# Make an outlet (Groundwater as a boundary condition)
# Add again later: groundwater=p.NewOutlet(name = 'outlet', x=0, y=0, z=-1) # creates a Dirichlet boundary condition with the potential z and add it to the list of outlets

# TODO make groundwater impermeable or remove it

# Simulation period
start = cmf.Time(1, 1, 2010)
end = start + timedelta(days=250)  # run cmf for 250 days
dt = cmf.h  # time step

# Make solver
solver = cmf.CVodeIntegrator(p, 1e-6)
solver.t = start

# Setup
initial_c_2H_soil = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
initial_c_18O_soil = iso.delta_to_concentration(delta_i=-8, solute_i="18O")
initial_c_2H_atmosphere = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
initial_c_18O_atmosphere = iso.delta_to_concentration(delta_i=-8, solute_i="18O")

T_atmosphere = 303.0
rH_atmosphere = 0.20

my_iso_atmosphere = iso_atmosphere(
    initial_c_atmosphere={"2H": initial_c_2H_atmosphere, "18O": initial_c_18O_atmosphere},
    # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
    initial_T_atmosphere=T_atmosphere,  # temperature in the atmosphere in Kelvin
    initial_Rh_atmosphere=rH_atmosphere)  # relative humidity of the atmosphere (-)

# my_iso_cell = iso_cell(atmosphere=my_iso_atmosphere, area=1.0)

i_upper_boundary = 0.0
my_layers = []
for i_lower_boundary, i_cmf_layer in zip(lower_boundaries_of_layer, c.layers):
    my_iso_layer = iso_layer(upper_boundary=i_upper_boundary,
                             lower_boundary=i_lower_boundary,
                             initial_c_solutes={"2H": initial_c_2H_soil, "18O": initial_c_18O_soil},
                             # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                             initial_theta=i_cmf_layer.theta,  # volumetric water content m3/m3
                             theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                             theta_sat=i_cmf_layer.porosity,  # volumetric water content m3/m3 at saturation
                             initial_T=303,  # soil temperature in Kelvin
                             porosity=i_cmf_layer.porosity,  # porosity of the soil m3/m3
                             tortuosity=0.67)  # tortuosity of the soil m/m
    i_upper_boundary = i_lower_boundary
    my_layers.append(my_iso_layer)

# The run time loop. Iterates over the outer timestep of the model
theta_layers = []
ql_up_t1 = []
ql_down_t1 = []
evap = []
for t in solver.run(solver.t, end, dt):

    # get fluxes of the next time step
    ql_layers_up_t1 = []  # flux in m3/d during the last time step
    ql_layers_down_t1 = []  # flux in m3/d during the last time step

    # Todo: Implement like this:
    for i_iso_layer, i_cmf_layer in zip(my_layers, c.layers):
        if i_cmf_layer.upper is None and i_cmf_layer.lower is not None:
            ql_layers_up_t1.append(0.0)
            ql_layers_down_t1.append(i_cmf_layer.flux_to(i_cmf_layer.lower, t))
        elif i_cmf_layer.upper is not None and i_cmf_layer.lower is not None:
            ql_layers_up_t1.append(i_cmf_layer.flux_to(i_cmf_layer.upper, t))
            ql_layers_down_t1.append(i_cmf_layer.flux_to(i_cmf_layer.lower, t))
        elif i_cmf_layer.upper is not None and i_cmf_layer.lower is None:
            ql_layers_up_t1.append(i_cmf_layer.flux_to(i_cmf_layer.upper, t))
            ql_layers_down_t1.append(0.0)

    ql_up_t1.append(ql_layers_up_t1)
    ql_down_t1.append(ql_layers_down_t1)
    theta_layers.append(c.layers.wetness * c.layers.porosity)

    # Evaporation
    E_pot = -1.005e-5  # kg/(m2*s)

    # calculate actual Evaporation for the next time step
    ha_surface = rH_atmosphere * iso_atmosphere.pv_sat(T_atmosphere) / iso_atmosphere.pv_sat(my_layers[0].T)

    actual_Evaporation = E_pot * (c.layers[0].wetness - ha_surface) / (
                1 - ha_surface) * c.area  # calculates the actual evaporation in kg/day
    E_act.flux = cmf.timeseries_from_scalar(actual_Evaporation)
    evap.append(actual_Evaporation)
    solver.reset()

theta_layers.append(theta_layers[-1])  # tehta value for last time step as the preceeding one

my_time = iso_time(final_time=250*24,
                   delta_time=1,
                   time_units='hours')  # 'seconds' , 'minutes', 'hours', 'days'

# Run isotope model
Ciso = {}
for solute in ['2H', '18O']:
    C_t = [[layer.c_solutes[solute] for layer in my_layers]]  # Initial isotopic concentration for all layers
    for t in range(my_time.time_steps-1):

        theta_t0 = theta_layers[t]
        theta_t1 = theta_layers[t+1]

        qv_up = [0] * len(my_layers)
        qv_down = [0] * len(my_layers)

        flux = [ql_up_t1[t], ql_down_t1[t], qv_up, qv_down, [theta_t0, theta_t1]]
        evaporation = evap[t]

        def BC(evap):
            E_pot = -1.005e-5  # kg/(m2*s)
            # Boundary conditions

            U_boundary_conc = iso.delta_to_concentration(delta_i=-0, solute_i='2H')
            L_boundary_conc = iso.delta_to_concentration(delta_i=-65, solute_i='2H')

            BC = BoundaryCondition()
            BC.upper_boundary('atmosphere', E_pot)  # neuman[flux] / dirichlet[constant conc] / atmosphere[potential evaporation]
            BC.lower_boundary('dirichlet', L_boundary_conc)  # TODO: need to check lower boundary

            return BC

        B_C = BC(evaporation)

        # Run simulations
        Ci_t = iso.run_1D_model(c, my_iso_atmosphere, my_layers, flux, B_C, solute,
                                ignore_alpha_i=ignorealphai,
                                ignore_alpha_i_k=ignorealphaik,
                                ignore_dl_i=ignoredl,
                                ignore_dv_i=ignoredv)

        updated_layers = iso.update_c_i(Ci_t, solute, my_layers)
        my_layers = updated_layers

        C_t.append(Ci_t)

    Ciso[solute] = C_t

C2H = np.array(Ciso['2H']).T
#C18O = np.array(Ciso['18O']).T

# Convert Result into delta notation
my_vectorized_function = np.vectorize(iso.concentration_to_delta, excluded=['solute_i'])
C_delta = my_vectorized_function(C2H, solute_i='2H')
#C18O_delta = my_vectorized_function(C18O, solute_i='18O')

#Visualize
plot = Visualize(my_layers, my_time)
plot.profile(C_delta, solute_i='2H', print_time_steps=24*20)
plot.breakthrough(C_delta, solute_i='2H', print_steps=2)

plt.plot(evap, label='evaporation')
plt.legend()
plt.show()

layer = 0
theta_l = [theta[layer] for theta in theta_layers]
plt.plot(theta_l, label='theta_top_layer')
plt.legend()
plt.show()

plt.plot(ql_up_t1[-1], range(len(ql_up_t1[-1])), label='flux at final time')
plt.legend()
plt.show()
