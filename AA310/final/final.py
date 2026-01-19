import sys 
sys.path.append('../')

import orbital_equations_of_motion as om 
import planetary_data as pd
import vector_functions as vf 

import math
import numpy as np
from matplotlib import pyplot as plt 
from matplotlib.patches import Circle

""" Question A """

print('\nA\n')

#Find the orbital parameters of the new glenn rocket booster

#Given 
alt = 252287 * 0.0003048  # Convert to ft to Km 
va = 0.00044704 * 4674    # Convert to mph to Km/s 
mu = pd.MU_EARTH_KM       # Mu Earth

ra = alt + pd.RADIUS_EARTH_KM 

v_circ = om.get_circular_orbit_speed(
  mu, ra
)

print(f'Vcirc at this alt {v_circ:.5g}, current speed {va:.5g}')

h = ra*va

#v_t = mu/h(1+ecos(theta)) Curtis 2.48
# At apoapsis theta = pi, cos theta = 1
e = 1 - (h*va/mu)

rp = ((1-e)/(1+e)) * ra  # Curtis 2.84
a = (rp + ra)/2          # Ellipitical Orbits Lecture

booster_state = om.orbital_state(
  a, e, mu, km=True 
)

om.print_state(
  booster_state, label='Booster'
)

""" Question B """

print('\nB')

# Sitting in the lagrange point 

dist_to_earth = pd.SUN_EARTH_KM
dist_to_l2 = dist_to_earth + (1.5 * (10**6))

fig, ax = plt.subplots()

ax.set_xlabel("Km")
ax.set_ylabel("Km")

ax.plot(0,0, 'yo', label='Sun')
ax.plot(
  dist_to_earth * math.cos(math.pi/4), dist_to_earth * math.sin(math.pi/4), 'bo', label='Earth')
ax.plot(
  dist_to_l2 * math.cos(math.pi/4), dist_to_l2 * math.sin(math.pi/4), 'go', label ='L2')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=3)

ax.set_aspect('equal')

# plt.show()

#Finding the eccentricity

rp = ra            #Booster ra
ra = 1.5 * (10**6) #km 

e = (ra - rp)/(ra+rp) # Curtis 2.84
a = (rp+ra)/2         # Semi-Major

mu = pd.MU_EARTH_KM

l2_state = om.orbital_state(
  a,e, mu, km=True
)

om.print_state(
  l2_state, label='Lagrange Orbit'
)

va # defined from earlier question

print(f'\nNeeded Dv {l2_state['v_p (km/s)'] - va:.5g}')

""" Question C """

# Assuming that we maintain the same speed at 200 km that the booster delievered
v_esc = om.get_v_esc(
  mu, pd.RADIUS_EARTH_KM + 200
)

# change in inclination 
i = 6.05 

print(f'\nEscape Velocity {v_esc:.5g}')

# Determine the plane change
dv_plane_change = om.inclination_change(
  v_esc, i
)

print(f'Plane Change Dv {dv_plane_change:.5g}')

va # From the first question 

dv_plane_change_stage_1 = om.inclination_change(
  va, i
)

print(f'Plane Change Dv {dv_plane_change_stage_1:.5g} km/s')

""" Question D """

# Get the transfer orbit from Earth to Mars

rp = pd.SUN_EARTH_KM + (1.5 * (10**6))
ra = pd.SUN_MARS_KM 

e = (ra - rp)/(rp + ra) # Curtis 2.84
a = (rp + ra)/2         # Ellicpitial Orbit Lecture Slides

xfer_state = om.orbital_state(
  a,e,pd.MU_SUN_KM, km=True
)

om.print_state(
  xfer_state
)

T_mars = om.get_period(pd.SUN_MARS_KM, pd.MU_SUN_KM)
w_mars = (2*math.pi)/T_mars

T_earth = 365 * 3600 * 24
w_earth = (2*math.pi)/T_earth

# Get angular distance covered by mars in one Earth year
# Add (2pi - this to the total distance, earth "gains" this much
# phase per year on mars)
phi_year = (math.pi * 2) - (w_mars * T_earth)

print(f'\nPhase Difference Created by waiting a year {np.rad2deg(phi_year):.5g} Deg')

# Get phase difference created while waiting for hohmann transfer
phi_hohmann = om.get_phase_difference_in_out(
  xfer_state['T (s)'], w_mars
)

print(f'\nPhase Difference Created by hohmann {np.rad2deg(phi_hohmann):.5g} Deg')

l2_state # from previous question

# Phase difference created during time getting to L2
phi_l2 = (w_earth * l2_state['T (s)']/2) - (w_mars * l2_state['T (s)']/2)

print(f'\nPhase Difference Created by waiting for L2 {np.rad2deg(phi_l2):.5g} Deg')

# Phase difference created during the 3 minutes and 10 seconds of stage 1
phi_launch = (w_earth * 190) - (w_mars * 190)

print(f'\nPhase Difference Created by waiting for launch {np.rad2deg(phi_launch):.5g} Deg')

phi = phi_launch + phi_l2 + phi_year + phi_hohmann

print(f'\nPhase Difference at time of New Glenn Launch {np.rad2deg(phi):.5g} Deg')

""" Question E """

delta = math.pi/2

# Turning angle is the change in angle between V_inf after a planetary fly by
e = om.get_eccentricity_turn_angle(delta, deg=False)
print(f'Eccentricity of the Hyperbolic Orbit {e:.5g}')

# Given 
rp = 40000 + pd.RADIUS_EARTH_KM # km

mu = pd.MU_EARTH_KM

flyby = om.hyperbolic_state(
  e, rp, mu
)

om.print_state(flyby, label='Hyperbolic')

# Get the transfer orbit from Earth to Mars
rp = pd.SUN_EARTH_KM
ra = pd.SUN_MARS_KM 

e = (ra - rp)/(rp + ra) # Curtis 2.84
a = (rp + ra)/2         # Ellicpitial Orbit Lecture Slides

xfer_state = om.orbital_state(
  a,e,pd.MU_SUN_KM, km=True
)

om.print_state(
  xfer_state
)

# Speed leaving Earth in heliocentric frame
v_sat_sun  = flyby['V_inf (km/s)'] + pd.V_EARTH_SUN_KM

# Speed required to enter mars transfer orbit
v_required = xfer_state['v_p (km/s)']

print(f'Speed leaving Earth {v_sat_sun:.5g} km/s, Speed required to enter Xfer orbit {v_required:.5g} km/s')

dv_req = v_required - v_sat_sun 

print(f'Required DV to enter Xfer orbit {dv_req:.5g} km/s')

