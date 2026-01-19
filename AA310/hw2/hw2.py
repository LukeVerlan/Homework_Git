import sys 
sys.path.append('../')

import math

import orbital_equations_of_motion 
import vector_functions


""" Part 1: CubeSat"""
r = [7000, 0, 0] #km
v = [0, 8, 0] #km/s
mass = 6 #kg

h = orbital_equations_of_motion.spec_ang_momentum_from_v_and_r(v,r)
print(f' Specific Angular Momentum Vector {h} km^2/s \n'
      f' Angular Momentum Vector {[direction * mass for direction in h]} kg * km^2 /s\n'
      f' Angular Momentum Magnitude {vector_functions.magnitude(h) * mass} kg * km^2 /s\n')


#Parts C D E 
r_vals = [[3500, 6062, 0], [0, 7000, 0], [3500, -3500, 0]]
v_vals = [[-6.928, 4.000, 0], [8, 0, 0], [8, 8, 0]]

for i in range(len(r_vals)):

  r = r_vals[i]
  v = v_vals[i]

  h = orbital_equations_of_motion.spec_ang_momentum_from_v_and_r(v,r)
  print(f' Specific Angular Momentum Vector {h} km^2/s \n'
        f' Angular Momentum Vector {[direction * mass for direction in h]} kg * km^2 /s\n'
        f' Angular Momentum Magnitude {vector_functions.magnitude(h) * mass} kg * km^2 /s\n')
  
print('Part 2:')

""" Part 2: Wobble """
# in km, relative to the given orbit
#values pulled from Curtis Appendix A 
wobbles = {
  'Sun to Earth'    : ('Sun', 'Earth', 149.6*math.pow(10,6)), #semi major axis
  'Sun to Jupiter'  : ('Sun', 'Jupiter', 778.8*math.pow(10,6)), #semi major axis
  'Earth to Moon'   : ('Earth', 'Moon', 384.4*math.pow(10,3)), #semi major axis
  'Earth to CubeSat': ('Earth', 'CubeSat', orbital_equations_of_motion.altitude_to_orbital_radius(450)),
  'Earth to ISS'    : ('Earth', 'ISS', orbital_equations_of_motion.altitude_to_orbital_radius(400)) 
}

#in kg
masses = {
  'Sun'     : 1.989*math.pow(10,30),
  'Earth'   : 5.974*math.pow(10,24),
  'Moon'    : 73.48*math.pow(10,21),
  'CubeSat' : 3,
  'ISS'     : 510000,
  'Jupiter' : 1.899*math.pow(10,27)
}

for wobble in wobbles:

  (body1, body2, radius) = wobbles[wobble]
  m1 = masses[body1]
  m2 = masses[body2]

  print(f'Wobble relative to larger mass: {wobble}. \n'
        f'Wobble abs dist (km) {orbital_equations_of_motion.wobble(m1, m2, radius)}\n'
        f'Wobble percent diff (%) {orbital_equations_of_motion.wobble(m1, m2, radius, percent_diff=True)}\n\n')
  

""" Part 3: Keplers equations """

theta_vals = [0, math.pi/2, -math.pi/2, math.pi]
e_vals = [0, 0.5, 3, 1]
h = 1 
mu = 1

for e in e_vals:
  print(f'For e of {e}: ')
  for theta in theta_vals:
    radius, velocity = orbital_equations_of_motion.keplers_scalar_state(h, mu, e, theta)
    print(f'Theta: {round(theta,4)} radius: {radius} m; velocity: {velocity}')
"""Part 4: The holy grail of orbital mechanics"""

e  = 0.91612
a  = 35674000 #m
mu = 4.90 * math.pow(10,12)

orbital_state = orbital_equations_of_motion.orbital_state(a,e,mu)
orbital_equations_of_motion.print_state(orbital_state)

RADIUS_MOON = 1737000 #m
r_p = 50000 + RADIUS_MOON
r_a = 3000000

e = (r_a-r_p)/(r_a+r_p)
a = r_p/(1-e)

print(f'\n {e} {a} {mu}')

orbital_state = orbital_equations_of_motion.orbital_state(a,e,mu)
orbital_equations_of_motion.print_state(orbital_state)








