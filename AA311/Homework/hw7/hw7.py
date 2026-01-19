import sys
sys.path.append('../../')
import tools

import math

""" Q1 """

print('Q1')

#given 
S = 47 * (3.28**2) #ft^2
V = 250 * 1.4667 #ft/s
rho = 0.002378 #slugs/ft3
q = tools.dynamic_free_stream(rho, V)
cl = 1.2
L = tools.lift_drag_formula(cl, q, S)
W = 103047 * 0.22481 #lbs
n = L/W
g = 32.2 #ft/s2

R = (V**2)/(g*math.sqrt((n**2)-1))
w = (g*math.sqrt((n**2)-1)/V)

print(f'Rmin (m): {R/3.28:.5g}, Wmax (rad/s): {w:.5g}')

print('Q2')

#given
S = 987 #ft2
W = 25000 #lb
AR = 9.14 
e = 0.7
n = 0.8
Vmax = tools.mph_to_fts(229) #ft/s
rho = 1.8975 * pow(10,-3) #slugs/ft3
Pmax = tools.horses_to_ft_lb_s(1200) #ft lb/s 

Pa = 2 * Pmax * n 
D = Pa/Vmax
q = tools.dynamic_free_stream(rho, Vmax) 
cl = W/(q*S)
cd = D/(q*S)
cdi = tools.induced_drag_coef(cl, e, AR)
cd0 = cd - cdi
print(f'Cd0: {cd0:.5g}')

""" Q3 """

print('Q3')

c_mcg = 0.0050
Cl = 0.50
static_margin = 0.03 

c_mac = c_mcg - Cl*(static_margin)

print(f'Cmac: {c_mac:.5f}')

""" Q4 """

print("Q4")

V = 100 #m/s
S = 1.5 #m2 
c = 0.45 #m
L = 3675 #N
Mcg0 = -12.4 #N*m
Mcg = 20.67 #N*m
rho = 1.225 #kg/m3

q = tools.dynamic_free_stream(rho, V)
Cmcg0 = tools.get_moment_coef(Mcg0, q, S, c)
Cmcg = tools.get_moment_coef(Mcg, q, S, c)

Cmac = Cmcg0 

cl = tools.get_lift_drag_coef(L, q, S)

sm = (Cmcg - Cmac)/cl
print(sm)