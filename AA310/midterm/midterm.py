import sys
sys.path.append('../')

import planetary_data
import orbital_equations_of_motion
import vector_functions

import math
from matplotlib import pyplot as plt

""" Q 2 """

#A 

alt_ISS = 450 
alt_cube = 450.004
radius_ISS = orbital_equations_of_motion.altitude_to_orbital_radius_earth(alt_ISS, km=True)
print(f'Cube Orbital Radius in Km {radius_ISS:.3f}')

ra_cube = radius_ISS
rp_cube = radius_ISS - 50 #given

e = (ra_cube-rp_cube)/(ra_cube+rp_cube) # Curts 2.83
a = rp_cube/(1-e) # Curtis 2.73

cubesat_state = orbital_equations_of_motion.orbital_state(a, e, planetary_data.MU_EARTH_KM, km=True)
orbital_equations_of_motion.print_state(cubesat_state, label='Cubesat Orbit')

v_ISS = math.sqrt(planetary_data.MU_EARTH_KM/radius_ISS) # Curis 2.61
delta_v = cubesat_state['v_a (km/s)'] - v_ISS #Final - initial
print(f'\nDelta V in km/s {delta_v:5g}')

mu = planetary_data.MU_EARTH_KM

h_ISS = v_ISS*radius_ISS
# ax = orbital_equations_of_motion.orbit_plot(h_ISS, mu, 0, color='b', label='ISS')
# orbital_equations_of_motion.orbit_plot(
#   cubesat_state['h (km^2/s)'], mu, e, ax=ax, color='g', label='Cubesat')

# plt.show()

""" Q E """

T_ISS = 2 * math.pi * math.sqrt((radius_ISS**3)/mu) #2.83
#one orbit is 2pi radians, in this equation that divides out to be this equation 
#since we are doing a "phasing maneuver" of one orbit in 18 ISS orbits
T_sat2 = T_ISS - (T_ISS/18) #Phasing maneuver equation
T_dif = T_sat2 - T_ISS 
print(f'\nSigned difference in periods is {T_dif:3g} s, or {(T_dif/60):5g} min')

mu = planetary_data.MU_EARTH_KM
a_sat2 = (((T_sat2 / (2 * math.pi)) ** 2) * mu) ** (1/3)

#Because this orbit shares an apoapsis with the ISS' orbit, know that r_a sat2 = r_iss
r_a_sat2 = radius_ISS
e_sat2 = (r_a_sat2/a_sat2) - 1 # ra = a(1+e)
sat2_state = orbital_equations_of_motion.orbital_state(a_sat2, e_sat2, mu, km=True)

delta_v_sat2 = sat2_state['v_a (km/s)'] - v_ISS
print(f'Signed Delta V for sat 2 to be in this orbit {delta_v_sat2:5g} (km/s)')

#Because all sats start at apoapsis, we need to add half of the respective satellite period on
#top of the period of the ISS to account for starting from apo
#finds a signed difference of state2 - state1
def find_difference_in_theta(state1, state2, given_time):
  state1_e = state1['e']
  state1_T = state1['T (s)']
  state2_e = state2['e']
  state2_T = state2['T (s)']

  #find true anom using newton rhapson                    account for starting at apo
  state1_theta = orbital_equations_of_motion.theta_from_t(
    given_time + (state1_T/2), state1_e, state1_T)
  state2_theta = orbital_equations_of_motion.theta_from_t(
    given_time + (state2_T/2), state2_e, state2_T)

  return state2_theta - state1_theta

mu = planetary_data.MU_EARTH_KM
ISS = orbital_equations_of_motion.orbital_state(radius_ISS, 0, mu, km=True)
sat1 = cubesat_state #defined in first part of question 2
sat2 = sat2_state #defined in previous question 
T_ISS #defined in previous questions 

combinations = [(('ISS','Sat1'),ISS, sat1),
                (('ISS','Sat2'), ISS, sat2),
                (('Sat1','Sat2'),sat1,sat2)]

for combination in combinations:
  names, state1, state2 = combination
  name1, name2 = names
  print(f'The difference in true anom from {name1} to {name2} is '
        f'{find_difference_in_theta(state1, state2, T_ISS):5g} deg')

""" Q 3 """
#given
mu = planetary_data.MU_EARTH_KM

#because we only know T, we can only find a and epsilon
T = 48 * 3600 #hours to seconds 
a = (((T / (2 * math.pi)) ** 2) * mu) ** (1/3) #2.83
epsilon = - mu/(2*a)

# print(f' a: {a:5g} km, epsilon: {epsilon:5g} MJ/kg')

#given
alt_apo = 90000 #km
ra = orbital_equations_of_motion.altitude_to_orbital_radius_earth(alt_apo, km=True)

e = ra/a - 1 # ra = a(1+e)

state = orbital_equations_of_motion.orbital_state(a, e, mu, km=True)

T = state['T (s)']
theta = 90
t_to_90 = orbital_equations_of_motion.t_from_theta(theta, T, e)
t_south_pole = 2 * t_to_90
t_north_pole = T - t_south_pole
print(f'Time spent around North pole in hours {(t_north_pole/3600):.5g}')
print(f'Time spent around South pole in hours {(t_south_pole/3600):.5g}')

T, e #from previous problem

range = 42

t_cone = orbital_equations_of_motion.time_spent_in_cone(range, e, T)
print(f'Time spent with sensors seeing the north pole {(t_cone/3600):.5g} hours')

mu = planetary_data.MU_EARTH_KM

#satellite state defined earlier 
state

#circular deploy orbit
alt_deploy = 38000 #km
r_deploy = orbital_equations_of_motion.altitude_to_orbital_radius_earth(alt_deploy, km=True)

#for the transfer starting at apoapsis
deploy_state = orbital_equations_of_motion.orbital_state(r_deploy, 0, mu, km=True)

delta_v_apo, transfer_state_apo = orbital_equations_of_motion.hohmann_transfer(
  deploy_state, state, mu, start_peri=False
)

#for the transfer starting at periapsis
delta_v_peri, transfer_state_peri = orbital_equations_of_motion.hohmann_transfer(
  deploy_state, state, mu, start_peri=True
)

orbital_equations_of_motion.print_state(transfer_state_apo, label='starting from apoapsis')
orbital_equations_of_motion.print_state(transfer_state_peri, label='starting from periapsis')

print(f'Total magnitude Delta V, transfer starting at Apo {delta_v_apo:.5g} km/s')
print(f'Total magnitude Delta V, transfer starting at Peri {delta_v_peri:.5g} km/s')

