import sys
sys.path.append('../')

import orbital_equations_of_motion
import planetary_data
import vector_functions

"""Q 2"""

#perifocal in km/s
vr1_a = [-3, 7]
vab_a = [-2, 10]

vab_b = [-6, -3]
vr2_b = [-4, -4]

#compute respective delta v
delta_transfer_a = orbital_equations_of_motion.perifocal_delta_v(vr1_a, vab_a)
delta_transfer_b = orbital_equations_of_motion.perifocal_delta_v(vab_b, vr2_b)
delta_transfer_total = delta_transfer_a + delta_transfer_b 

print(f'Total delta transfer in km/s {delta_transfer_total:.4f}')

"""  B  """

import math

mu = planetary_data.MU_EARTH_KM
vp_a = vector_functions.magnitude(vr1_a) #Velocity at A
rp_ab = (mu)/(vp_a**2) # v_t = sqrt(mu/r_c) (2.63)
vp_ab = vector_functions.magnitude(vab_a)
h = vp_ab*rp_ab #(2.31)
e = (h**2)/(mu*rp_ab) - 1 #keplers 2nd law 
#given theta = 120 deg
theta = 120 * math.pi/180 
E = 2*math.atan(math.sqrt((1-e)/(1+e))*math.tan(theta/2)) # (3.13b)
Me = E - e * math.sin(E) #(3.14)
a = rp_ab/(1-e) #(2.73)
T = 2*math.pi*math.sqrt((a**3)/mu) # (2.83)
t = (Me/(2*math.pi))*T/60 # (3.15) time in minutes

print(f'Time to point B in minutes {t:.3f}')
 
""" C """
#circular orbit
time = 10*3600+(t*60) #hours to seconds maintaining the previous period
# theta = orbital_equations_of_motion.theta_from_t(time, e, T)
print(f'True anomaly after 10 hours of idle: {theta:.2f} deg')

""" 3 A """

print('\n A \n')

alt_ISS = 450000 #m
mu = planetary_data.MU_EARTH
r = orbital_equations_of_motion.altitude_to_orbital_radius_earth(alt_ISS,km=False)
T = 2 * math.pi * math.sqrt((r**3)/mu)
#given
theta = -14 #less than 20 degrees is valid for the small sine approximation
orbital_equations_of_motion.phasing_maneuver(theta, r, mu, T, label='Part A')


""" B """

print('\n B \n')

theta = 14
#given
t = 2 * 3600 * 24 #hours to seconds 
orbital_equations_of_motion.phasing_maneuver(theta, r, planetary_data.MU_EARTH, t, label='Part C')


""" 4 """

""" A """

from matplotlib import pyplot as plt 

#given
rp_alt = 1650 #km
rp = orbital_equations_of_motion.altitude_to_orbital_radius_earth(rp_alt, km=True)
e  = 0.58

a = rp/(1-e) #(2.73)

mu = planetary_data.MU_EARTH
r_meo = orbital_equations_of_motion.altitude_to_orbital_radius_earth(22000, km=True)
failState = orbital_equations_of_motion.orbital_state(
  a*1000, e, mu) #needs SI units

mu = planetary_data.MU_EARTH_KM
v_meo = math.sqrt(mu/r_meo) #(2.61)

h_meo = v_meo * r_meo #km^2/s

h_fail = failState['h (m^2/s)'] / (10**6) #put in km

#get ax, fail orbit in blue
ax = orbital_equations_of_motion.orbit_plot(h_fail, mu, e, label='Fail')

#meo orbit in red
orbital_equations_of_motion.orbit_plot(h_meo, mu, 0, ax, color='r', label='MEO')

# plt.show()

""" B """

print('\n 4 \n')

print('A \n')
#Need to find theta of transfer happens when r_meo = r_fail

#e, mu (km), r_meo (km) from previous problem
theta = math.acos((1/e)*((h_fail**2)/(mu*r_meo)-1))
print(f'Theta {theta * (180/math.pi):.2f} deg')

vr_fail = (mu/h_fail)*e*math.sin(theta) #(2.49)
vperp_fail = (mu/h_fail)*(1+e*math.cos(theta)) #(2.48)

#v_meo found in previous question

#vectors in vr, v_perp
vector_meo = [0, v_meo]
vector_fail = [vr_fail, vperp_fail]

# print(vr_fail)
# print(vperp_fail)
# print(vector_functions.magnitude(vector_fail))
print(v_meo)
delta_v = vector_functions.subtraction(vector_fail,vector_meo)
print(f'Delta V vector in format of Vr, Vperp  {[f"{v:.3f}" for v in delta_v]}, \n' 
      f'Delta V magnitude (km/s) {vector_functions.magnitude(delta_v):.3f}')

""" C """

print('\n c \n')

gamma_fail = orbital_equations_of_motion.flight_path_angle(
  vector_fail[1], vector_fail[0])
gamma_meo  = orbital_equations_of_motion.flight_path_angle(
  vector_meo[1], vector_meo[0])
delta_gamma = gamma_meo - gamma_fail

print(f'Delta gamma {delta_gamma:.3f} deg')

""" E """

#from question A we have fail state and r_meo already
failState, r_meo

mu = planetary_data.MU_EARTH #si units

meo_state = orbital_equations_of_motion.orbital_state(
  r_meo*1000, 0, mu) #correct MEO to m

orbital_equations_of_motion.print_state(failState)
orbital_equations_of_motion.print_state(meo_state)

delta_v_start_peri = orbital_equations_of_motion.hohmann_transfer(
  failState,meo_state,mu,km_s=True, start_peri=True)
delta_v_start_apo  = orbital_equations_of_motion.hohmann_transfer(
  failState,meo_state,mu,km_s=True, start_peri=False)

print(f'Hohmann Delta V starting from Peri {delta_v_start_peri:.4f} km/s\n'
      f'Hohmann Delta V starting from Apo  {delta_v_start_apo:.4f} km/s')


