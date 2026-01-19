import math 

#returns y when given y_1, y_2, x_1, x_2, and x 
# takes form y = ((y2-y1)/(x2-x1)) * (x-x1)
def interpolation(y1, y2, x1, x2, x):
  return (y1 + (((y2-y1)/(x2-x1))*(x-x1)))

#expects free stream values
def dynamic_free_stream(rho, velo):
  return (1/2) * rho * (velo**2)

#expects free stream values
def induced_drag_coef(cl, e, AR):
  return (cl**2)/(math.pi * e * AR)

def get_moment_coef(M, q, S, c):
  return M/(q*S*c)

def get_lift_drag_coef(L_or_D, q ,S):
  return L_or_D / (q * S)

#expects free stream
def lift_drag_formula(cl_or_cd, q, S):
  return cl_or_cd * q * S

def mph_to_fts(V):
  return V * 1.4667

def horses_to_ft_lb_s(P):
  return P * 550

def deg2rad(deg):
  return math.pi * deg/180 

def rad2deg(rad):
  return rad *180/math.pi