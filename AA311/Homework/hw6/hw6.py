import sys
sys.path.append('../../')

import tools
import math

""" 6.1 """

W = 5000
L = W
rho = 0.002377
v = 200  * (5280/3600) # ft/s
S = 200 
AR = 8.5
e = 0.93

q = tools.dynamic_free_stream(rho, v)

cl = L/(q*S)

cd = 2 * tools.induced_drag_coef(cl, e, AR)

D = L*cd/cl

""" 6.9 """

L_D = 7.7

theta = math.atan(1/7.7) 

print(f'Theta = {theta  * 180/math.pi}')

delta_x = 5000/math.tan(theta)

print(f'Delta X = {delta_x:.5g} ft')

""" 6.11 """

cd0 = 0.025 
AR = 6.72 
e = 0.9 

cl = math.sqrt(cd0 * math.pi * e * AR)

cd = 2 * cd0

L_D_max = cl/cd

print(f'L/D max = {L_D_max:.5g}')

""" Problem 6.16 """

W = 103047 #N
e = 0.87

AR = 6.5
S = 47
cd0 = 0.032
T = 2 * 40298
h = 5 / 3.28084 
Cl_max = 0.8
rho = 1.225 
mu = 0.02

b = math.sqrt(AR * S)
phi = ((16*h/b)**2)/(1 + ((16*h/b)**2))

cd = cd0 + (0.8**2) * phi / (math.pi * e * AR)

print(h)
print(phi)
print(cd)

vstall = math.sqrt(2 * W/(Cl_max * S * rho))

v_avg = vstall* 1.2 * 0.7
print(v_avg)

Davg = cd * (1/2) * (v_avg**2) * rho * S 
Lavg = Cl_max * (1/2) * (v_avg**2) * rho * S

Slo = (1.44 * W**2)/(9.81 * S * Cl_max * rho * (T - (Davg + mu * (W - Lavg))))
print(Slo)