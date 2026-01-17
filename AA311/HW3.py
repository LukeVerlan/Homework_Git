import math

R_cr = 5 * (10**5)
mu = 1.789 * pow(10, -5)
rho = 1.225 
L = 4 

def drag(v):
    x_cr = (R_cr * mu)/(rho * v)
    q = (rho * (v**2))/2
    R_l = (rho * v * L)/mu 
    c_len_turb = 0.074/(R_l**(1/5))
    drag_turb_len = q * c_len_turb * L * L
    c_a_T = 0.074/(R_cr**(1/5))
    drag_turb_a = q * L * x_cr * c_a_T
    c_a_L = 1.328/(R_cr**(1/2))
    d_b = drag_turb_len - drag_turb_a
    d_a = L * q * (x_cr * c_a_L)
    return 2 * (d_b + d_a)

n = (math.log(drag(40)/drag(20))/math.log(40/20))

gamma = 1.4
T = 248.55
rho = .6523
mu = 1.52 * (10**(-5))
R = 287

L = 2.89

a = math.sqrt(gamma * T * R)
M = 0.6

v = M*a

r_l = (v*rho*L)/mu

delta = 0.37 * L/(r_l**(1/5))

print(delta)


