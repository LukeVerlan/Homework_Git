import sys
sys.path.append('../')

import orbital_equations_of_motion
import planetary_data
import vector_functions

""" Q2 """

omega = 90
i = 90 
w = -45
rp = 12000
e = 0.5
theta = 45 

a = rp/(1-e)

orbital_equations_of_motion.geocentric_3d_plot(omega, i, w, theta, a, e)

""" Q3 """

r = [1000, 2000, 20000]
v = [3, 3, -0.3]

geo_state = orbital_equations_of_motion.get_geocentric_orbital_params(
  r, v
)

orbital_equations_of_motion.print_state(geo_state, label='Geocentric')

""" Q4 """

#given 
A = 8 
Ar = 0.05
a = 24
ar = 0.02
p = 8250 
pr = -0.25
H = 220/3.28
sidr = 36
lat = 47.660503

r,v = orbital_equations_of_motion.alg_5_4(
  A, Ar, a, ar, p, pr, lat, sidr, H
)

vec_units = ['I','J','K']


vector_functions.print_vector_formatted(r, vec_units, units='km', label='r')
vector_functions.print_vector_formatted(v, vec_units, units='km/s',label='v')

geo_state = orbital_equations_of_motion.get_geocentric_orbital_params(
  r, v
)

orbital_equations_of_motion.print_state(geo_state, 'Geocentric')




