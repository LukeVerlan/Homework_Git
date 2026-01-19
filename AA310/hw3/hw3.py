import sys
sys.path.append("../")

import planetary_data
import vector_functions
import orbital_equations_of_motion
import math 
from matplotlib import pyplot as plt
import numpy as np 


""" PART 1 """

# a 
def part_a(n):
  return math.sqrt((1-(1/n**2)))

n_vals = [1,2,3,4,5]

for n in n_vals:
  print(f'For n of {n}: e = {part_a(n)}')

"""Part b"""
a = 1
theta = np.linspace(0, 2 * np.pi, 500)

eccentricities = [part_a(n) for n in n_vals]

fig, ax = plt.subplots()
for e in eccentricities:
    r = a * (1 - e**2) / (1 + e*np.cos(theta))
    x = r * np.cos(theta) #Add + a*e for geometric center, otherwise plts focal pt
    y = r * np.sin(theta)
    ax.plot(x, y, label=f"e={e:.3f}")

ax.set_aspect('equal')
ax.set_title("T = const, varying eccentricity")
ax.legend(loc='lower right')
plt.show()

"""Part D"""

# #problem statement
# r = 1
# mu = 1
# v_1 = math.sqrt(mu/r)

# # returns a velocity ratio v_p/v_1 as a function of e 
# def velo_ratio_peri(e):
#     return math.sqrt((1+e)/(1-e)) #resolved vis-viva of vp/v1

# # returns a velocity ratio v_a/v_1 as a function of e 
# def velo_ratio_apo(e):
#     return math.sqrt((1-e)/(1+e)) #resolved vis-viva of va/v1

# #init lists
# vp_over_v1 = []
# va_over_v1 = []
# e_vals = []

# fig, ax = plt.subplots()

# for e in np.arange(0, 0.90001, .125):
#     e_vals.append(e)
#     vp_over_v1.append(velo_ratio_peri(e))
#     va_over_v1.append(velo_ratio_apo(e))

# plt.plot(e_vals, vp_over_v1, 'bo', label='Vp/V1')
# plt.plot(e_vals, va_over_v1, 'ro', label='Va/V1')
# plt.legend(loc='upper left')
# plt.xlabel('Eccentricity e')
# plt.ylabel('Velocity ratio')
# plt.title('Vp/V1 & Va/V1 as a function of e')
# plt.grid(True)
# plt.show()

""" Problem 2 """

# #returns the needed change in velocity to escape orbit 
# #takes in r as an orbital radius
# def needed_change_in_velo(r):
#   return math.sqrt(mu/r)*(math.sqrt(2)-1)

# orbital_radi ={
#   'LEO' : orbital_equations_of_motion.altitude_to_orbital_radius_earth(450 * (10**3)), #m from earth CoM
#   'MEO' : orbital_equations_of_motion.altitude_to_orbital_radius_earth(22200 * (10**3)), 
#   'GEO' : orbital_equations_of_motion.altitude_to_orbital_radius_earth(35700 *(10**3))
#   }

# for orbit in orbital_radi:
#   print(f'{orbit}: Delta V : {needed_change_in_velo(orbital_radi[orbit])} m/s')

# c = 2.998 * (10**8)

# r = (2*mu)/(c**2)
# mu = planetary_data.MU_EARTH

# v_new = math.sqrt(mu/orbital_radi['LEO']) + math.sqrt(mu/orbital_radi['LEO']) * (math.sqrt(2) - 1) / 2
# r = orbital_radi['LEO'] 
# h = v_new * r
# e = ((h**2)/(mu*r)) - 1
# a = r/(1-e)

# state = orbital_equations_of_motion.orbital_state(a, e, mu)

# orbital_equations_of_motion.print_state(state)

""" Problem 3 """

mu = planetary_data.MU_EARTH

rp = orbital_equations_of_motion.altitude_to_orbital_radius_earth(428000)
vp = 20880  # m/s

# h = rp * vp
# e = (h**2 / (rp * mu)) - 1

# state = orbital_equations_of_motion.hyperbolic_state(e, rp, mu)

# theta_inf = state['theta inf (deg)'] * math.pi / 180
# h = state['h (m^2/s)']
# e = state['e']
# a = state['a (m)']  # negative for hyperbola

# deg = 90
# theta = np.linspace(-(deg * math.pi / 180) + 0.01, (deg * math.pi / 180) - 0.01, 1000)

# def keplers_second_law(h, e, mu, theta):
#     return (h**2 / mu) / (1 + e * np.cos(theta))

# radius = keplers_second_law(h, e, mu, theta)

# # Shift the orbit so that the focus is at x = 0
# x_shift = a + rp
# x = radius * np.cos(theta) - x_shift
# y = radius * np.sin(theta)

# fig, ax = plt.subplots()

# # Center the plot around 0 with equal aspect ratio
# ax.set_aspect('equal', adjustable='datalim')

# # Automatically get the extents of the hyperbola
# x_min, x_max = np.min(x), np.max(x)
# y_min, y_max = np.min(y), np.max(y)

# # Make limits symmetric around zero
# x_range = max(abs(x_min), abs(x_max))
# y_range = max(abs(y_min), abs(y_max))

# ax.set_xlim(-x_range, x_range)
# ax.set_ylim(-y_range, y_range)

# # Add a grid for clarity
# ax.grid(True)

# ax.plot(x, y, label='Hyperbolic trajectory')
# ax.plot(0 - a - rp, 0, 'o', label='Focus (Earth)')
# ax.set_xlabel('x (m)')
# ax.set_ylabel('y (m)')
# ax.legend()
# ax.set_title(f'Hyperbolic Orbit (e={round(e,3)}, rp={round(rp,2)})')
# plt.show()

# e = 1.346 

# state = orbital_equations_of_motion.hyperbolic_state(e, rp, mu)
# orbital_equations_of_motion.print_state(state)
# print(state['V_p (m/s)'])

""" Part 4 """

angles = [0, 90, 135, -120] #deg
colors = ['r', 'b', 'g', 'y']

a = planetary_data.RADIUS_EARTH #m

e = 0.5 #statement

fig, ax = plt.subplots()

for i, angle in enumerate(angles):
  orbital_equations_of_motion.perifocal_plot(a,e,angle,ax, color=colors[i])

plt.show()

