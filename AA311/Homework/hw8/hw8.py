import sys 
sys.path.append('../../')

import tools

""" Q1 """

e0        = 0 
a_e       = 0.42
alpha_wb  = 5 #deg
it        = 2 #deg
St        = 0.4
c         = 0.45
S         = 1.5
l         = 1
at       = 0.12
v         = 100
rho       = 1.225
L         = 4134
cm0       = -0.003 
sm        = 0.02

q = tools.dynamic_free_stream(rho, v)

cl = L/(q*S)

a = cl/alpha_wb

e = e0 + a_e*alpha_wb

alpha_t = alpha_wb - it - e 

clt = (at)*(alpha_t)

a = cl/alpha_wb

Vh = (l*St)/(c*S)

cmcg = cm0 + cl*(0.02 - (Vh * (at/a)*(1-a_e)))+ (Vh*at*(it + e0))


Mcg = cmcg * q * S * c 

print("Total cm_cg =", cmcg)
print("Final moment about CG Mcg =", Mcg, "N·m")

""" Q2 """

h = 0.26
hac = 0.24 

hn = 0.24 + 0.59*(0.12/0.09)*(1-0.42)

print(hn-h)

""" Q3 """

F = 1 - (1/0.12)*(0.04*(-0.007/-0.012))

hpn = hac + F * Vh*(at/a)*(1-a_e)
print(hpn)

print((hpn-h)/(hn-h))

