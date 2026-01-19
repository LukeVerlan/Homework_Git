import sys
sys.path.append('../')

import planetary_data
import orbital_equations_of_motion
import vector_functions
import numpy as np 
from matplotlib import pyplot as plt

import math

""" Q 1 """

print('\n Question 1 \n')

#definitely dont wanna park here but it is the technical minimum
min_parking_alt = 0 

#max is at the edge of the soi (convert to an alt)
max_parking_alt = planetary_data.SOI_MARS_KM - planetary_data.RADIUS_MARS_KM

# alts = np.linspace(min_parking_alt, max_parking_alt, 10000)
# delta_vs = orbital_equations_of_motion.arrival_delta_v_mars(alts)

def delta_v_plot(alts, delta_vs):

  fig, ax = plt.subplots()
  ax.set_title('Alts Vs Delta Vs')
  ax.set_ylabel(f'Delta Vs (Km/s)')
  ax.set_xlabel(f'Alts (Km)')

  ax.semilogx(alts, delta_vs, color='g', label='Delta V vs Parking Alt')
  plt.show()

# delta_v_plot(alts, delta_vs)

# alts = np.linspace(min_parking_alt, max_parking_alt, 10000)
# delta_vs = orbital_equations_of_motion.three_delta_V(alts)

# delta_v_plot(alts, delta_vs)

# alts = np.linspace(min_parking_alt, max_parking_alt, 10000)
# delta_vs = orbital_equations_of_motion.fourth_delta_v(alts)

# delta_v_plot(alts, delta_vs)

# small_delta_v = orbital_equations_of_motion.fourth_delta_v(0)
# print(f'Smallest Delta V = {small_delta_v:.5g} Km/s')


aiming_radius_0 = orbital_equations_of_motion.aiming_radius_from_parking_radius_mars(
  0
)

aiming_radius_150 = orbital_equations_of_motion.aiming_radius_from_parking_radius_mars(
  150
)

print(f'Aiming Radius at 150 km {aiming_radius_150:.5g} Km')
print(f'Aiming Radius at 0 km {aiming_radius_0:.5g} Km')

""" Question 2 """

print("\n Question 2 \n")

planet_speed = planetary_data.V_MARS_SUN_KM
r_earth_sun = planetary_data.SUN_EARTH_KM
r_mars_sun = planetary_data.SUN_MARS_KM

v_inf = orbital_equations_of_motion.get_v_inf_arrival_outer(
  r_earth_sun, r_mars_sun, planet_speed)

print(f'V_inf entering mars {v_inf:.5g} Km/s')

#given
delta = 5000 #km
mu = planetary_data.MU_MARS_KM

a = mu/(v_inf**2) # Curtis 2.112
e = math.sqrt(((delta/a)**2) + 1) # Curtis 2.107

rp = a * (e-1) # curtis 2.101

print(f'Rp of the current hyperbolic orbit {rp:.5g} Km')

#part c 

#given 
r_target = 3800
a, e, mu #the same since last question
hyper_state = orbital_equations_of_motion.hyperbolic_state(e, rp, mu)
h = hyper_state['h (km^2/s)']

# Keplers 2nd r = h^2/mu(1/(1+ecos(theta)))
# Need to add negative because the angle is opposite of our sign convention
# for a hyperbola 
theta = -math.acos((1/e)*((h**2)/(r_target*mu) - 1))
print(f'Theta at which r_target is reached {theta * 180/math.pi:.5g} Deg')

#part d 

# Curtis 2.52 tan gamma = (esin(theta))/(1+ecos(theta))
gamma = math.atan(e*math.sin(theta)/(1 + (e * math.cos(theta)))) * 180/math.pi

print(f'Flight path angle at theta = {theta * (180/math.pi):.5g} deg, gamma = {gamma:.5g} deg')


mu, h, theta, e #from previous

# Curtis 2.49
v_r = (mu/h)*e*math.sin(theta)

# Curtis 2.48
v_perp = (mu/h)*(1 + (e*math.cos(theta)))

#want speed, not vector
v = math.sqrt((v_r**2) + (v_perp**2))

print(f'Speed at this point, |V| = {v:.5g} km/s')

#speed of the parking orbit
v_park = [math.sqrt(mu/r_target), 0]
v_hyper = [v_perp, v_r]

#get delta V vector
delta_v_vector = vector_functions.subtraction(v_park, v_hyper)

#find magnitude of that vector
delta_v = vector_functions.magnitude(delta_v_vector)

print(f'Magnitude of Delta V to change orbits {delta_v:.5g} km/s')

#Because these velocities were taken with respect to mars, 
#Can just compare the magnitude of these two vectors to find

mag_v_hyper = vector_functions.magnitude(v_hyper)

mag_v_park = vector_functions.magnitude(v_park)

#if |v_hyper| > |v_park| the satellite lost speed with respect to mars

if (mag_v_hyper > mag_v_park):
  print(f'Satelite lost {abs(mag_v_hyper - mag_v_park):.5g} km/s relative to Mars')
else:
  print(f'Satelite gained {abs(mag_v_hyper - mag_v_park):.5g} km/s relative to Mars')

# We that v_r makes a |theta| (found above) degree angle with mars orbit 
# relative to the sun. Therefore these vectors are 90 - |theta| degrees 
# rotated from what they should be relative to the sun. Using a rotation 
# matrix we can find the new magnitudes of these vectors after being rotated.
# Additionally, we need to add the velocity of mars to the perpendicular 
# components after rotating

def rotation_matrix(v_r, v_perp, theta):
  return [
    v_perp * math.cos(theta) - v_r * math.sin(theta),
    v_perp * math.sin(theta) + v_r * math.cos(theta)
  ]

#rotates below the horizontal
v_park_sun = rotation_matrix(v_park[1], v_park[0], -(math.pi/2 - abs(theta))) 
v_park_sun[0] += planetary_data.V_MARS_SUN_KM
v_hyper_sun = rotation_matrix(v_hyper[1], v_hyper[0], -(math.pi/2 - abs(theta)))
v_hyper_sun[0] += planetary_data.V_MARS_SUN_KM

mag_v_hyper_sun = vector_functions.magnitude(v_hyper_sun)
mag_v_park_sun = vector_functions.magnitude(v_park_sun)

if (mag_v_hyper_sun > mag_v_park_sun):
  print(f'Satelite lost {abs(mag_v_hyper_sun - mag_v_park_sun):.5g} km/s relative to Sun')
else:
  print(f'Satelite gained {abs(mag_v_hyper_sun - mag_v_park_sun):.5g} km/s relative to Sun')









""" Question 3 """

print('\n Question 3 \n')

r_mars = planetary_data.SUN_MARS_KM
r_earth = planetary_data.SUN_EARTH_KM
mu = planetary_data.MU_SUN_KM

state_mars = orbital_equations_of_motion.orbital_state(
  r_mars,0,mu,km=True)

state_earth = orbital_equations_of_motion.orbital_state(
  r_earth,0,mu,km=True
)

delta_v, transfer_state = orbital_equations_of_motion.hohmann_transfer(
  state_earth, state_mars, mu, start_peri=True
)

Txfer = transfer_state['T (s)']
t_apo = Txfer/2

print(f'Time to Mars from Earth: {t_apo / (3600 * 24):.5g} days')

#B 

r_mars, r_earth, mu, t_apo # same from previous

T_mars = 2 * math.pi * math.sqrt((r_mars**3)/mu)

n_mars = 2 * math.pi / T_mars

phi_i = orbital_equations_of_motion.hohmann_inner_outer_phase_difference(
  n_mars, t_apo
)

print(f'Initial Phase Difference Phi: {phi_i:.5g} deg')

Txfer, r_earth, mu #same from previous

t_peri = Txfer/2 

T_earth = 2 * math.pi * math.sqrt((r_earth**3)/mu)

n_earth = 2 * math.pi / T_earth

phi_f = orbital_equations_of_motion.hohmann_inner_outer_phase_difference(
  n_earth, t_peri
)

print(f'Return Phase Difference Phi: {phi_f:.5g} deg')

delta_phi = orbital_equations_of_motion.compute_delta_phi(
  n_earth, Txfer
)

t_wait = orbital_equations_of_motion.waiting_time_after_inner_outer(
  n_mars, n_earth, delta_phi
)

print(f'Time spent waiting {t_wait/(3600*24):.5g} days')

Total_time_for_hohmann = (Txfer + t_wait)/(3600 * 24)

print(f'Total time for ideal hohmann transfer {Total_time_for_hohmann:.5g} days')

""" Question 4 """

#c 

r_sun_earth = planetary_data.SUN_EARTH_KM
r_sun_jupiter = planetary_data.SUN_JUPITER_KM
speed_jupiter = planetary_data.V_JUPITER_SUN_KM

v_inf = orbital_equations_of_motion.get_v_inf_arrival_outer(
  r_sun_earth, r_sun_jupiter, speed_jupiter
  )

print(f'|V_inf| entering jupiter: {v_inf:.5g} km/s')

v_inf # same as previous 

# know |v_inf_in| = |v_inf_out| 
mu = planetary_data.MU_JUPITER_KM
a_hyper = mu/(v_inf**2) # Curtis 2.112

#Found the same way as previous velocities
V_jupiter = planetary_data.V_JUPITER_SUN_KM 
print(V_jupiter)
# V = math.sqrt(mu_sun/r_jupiter_sun)

#given exit properties 
V_sat_saturn = 14.410
gamma_saturn = -23.044 * (math.pi/180)

# Law of sines
turn_angle = math.asin((V_sat_saturn*math.sin(gamma_saturn)/v_inf))

print(turn_angle)

# Adjust for obtuse 
turn_angle = math.pi/2 + (math.pi/2 - abs(turn_angle))

print(f'Turn Angle: {turn_angle*(180/math.pi):.5g} deg') 

# Curtis 2.100 
e = 1/(math.sin(turn_angle/2))

# Curtis 2.107
delta = a_hyper * math.sqrt((e**2) - 1)

print(f'Required Aiming Radius {delta:.5g} km')




