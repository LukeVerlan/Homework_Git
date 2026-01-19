import sys
import math
sys.path.append('../')

import vector_functions
import orbital_equations_of_motion

A = [6,18,5]
B = [-2,-5,-8]
C = [2,1,12]

RADIUS_EARTH = 6371000 # m

# Question 1

"""
Use your functions to determine how aligned is vector A with the plane created
by B and C, comment on the result

A scalar triple product first creates a plane between the vectors in the second half of the
arguments, and then performs a dot product on that vector. Because the cross product returns an orthogonal
vector to the vectors in the cross product, the closer the absolute value of scalar triple product of a,b,c is to 0, 
the more allinged vector a is with the plane created by vector b & c 

"""

print(f"A * (B X C) {abs(vector_functions.scalar_triple_product(A,B,C))}")
print(f"C * (A X B) {abs(vector_functions.scalar_triple_product(C,A,B))}")
print(f"B * (A X C) {abs(vector_functions.scalar_triple_product(B,A,C))}")

print(vector_functions.magnitude(A) * 
      vector_functions.magnitude(vector_functions.cross_product(B,C)))

"""
When a vector is colinear with another vector the dot product will return a product of their 
magnitudes, ||A|| * ||B X C|| is ~ 1044, the scalar triple product returns 128 meaning that relative 
to a complete colinear vector, vector A is relatively close to the plane created by vectors B and C.
Additionally the scalar triple product is cyclic, meaning that the vectors can be in any orientation
and the value will be the same. Qualitatively the vectors are always the same distance apart and therefore
it doesnt matter what vectors you choose to be the plane. 

"""

leo_r = 460000 + RADIUS_EARTH # m
meo_r = 20200000 + RADIUS_EARTH # m 
geo_r = 35780000 + RADIUS_EARTH # m

print(f' LEO : {round(orbital_equations_of_motion.circular_orbit_velocity(leo_r) / 1000, 3)} km/s')
print(f' MEO : {round(orbital_equations_of_motion.circular_orbit_velocity(meo_r) / 1000, 3)} km/s')
print(f' GEO : {round(orbital_equations_of_motion.circular_orbit_velocity(geo_r) / 1000, 3)} km/s')

#Question 4 

# (v_0^2 + 2gr)/3gr = cos(phi)

v_0 = 0
lhs = 2/3
phi = math.acos(lhs) * 180/math.pi #compute angle
print(f'Departure Angle : {round(phi,3)}')

#v_0 = sqrt(gr(3cos(phi) - 2))

g = 9.81 # m/s 
phi = 45 * math.pi/180 #convert to radians
const = math.sqrt(g * (3*math.cos(phi) - 2)) #compute the solved version
print(f'V_0 : {round(const,3)} * r^(1/2)')