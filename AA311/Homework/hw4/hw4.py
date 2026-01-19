import sys
sys.path.append('../../')

import tools

""" Question 1 """

#imperial

rho = 2.3769 * pow(10,-3)
v = 100 
c = 3

c_l = 0.7
c_d = 0.0072
c_m = -0.025

q = tools.dynamic_free_stream(rho, v)

L = q * c * c_l
D = q * c * c_d
M = q * (c**2) * c_m

# print(q)

# print(f'Lift {L:.2f} Drag {D:.2f} Moment {M:.2f}')

""" Question 2 """

import math

c_l = 0.85
M_inf = 0.7

c_l_inc = c_l * math.sqrt(1-(M_inf**2))

# print(c_l_inc)

""" Question 3 """

T = 390.53 #R
R = 1716
Gamma = 1.4 
M = 2.2
rho = 7.1028 * pow(10,-4)
w = 36000
s = 210
a = math.sqrt(Gamma*R*T)
v = M*a

q = (1/2)*(rho)*(v**2)

c_l = w/(q*s)
alpha_rad = c_l * math.sqrt(M**2 - 1) / 4 
alpha_deg = alpha_rad * 180/math.pi
# print(alpha_deg)

c_dw = (4*(alpha_rad**2))/math.sqrt((M**2) - 1)

D_w = 2 * c_dw * q * s
print(c_dw)

""" Problem 6.16 """

W = 103047 #N
e = 0.87

AR = 6.5
S = 47
cd0 = 0.032
T = 2 * 40298
h = 5 

b = math.sqrt(AR * S)
phi = ((16*h/b)**2)/(1 + (16*h/b)**2)

print(phi)