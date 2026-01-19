import math 

a0 = 0.10625
e=0.9
ar = 5

a = a0/(1+((57.3*a0)/(math.pi*e*ar)))
alpha0 = -1.75
alpha = 6
cl = a*(alpha-alpha0)
# print(f'Cl = {cl:5g}')

Cp = 0.004
Cd = Cp + (cl**2)/(math.pi * e * ar)

# print(f'Cd = {Cd:5g}')

""" 5.28 """

rho = 0.002378
v = 21.8
q=(1/2)*rho*(v**2)
L = (1/16)
S = 1 
alpha = 3 #deg
c = 2 *  math.pi * alpha * (math.pi/180)

print(f'Q = {q:.5g}')

cl = L/(q*S)

print(f'Cl = {cl:.5g}')

cl_m = L/(q*c)

print(f'Cl_m = {cl_m:.5g}')
