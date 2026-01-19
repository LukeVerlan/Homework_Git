import math
import vector_functions
import planetary_data
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.mplot3d import Axes3D

"""
This file contains various orbital mechanics equations of motion 

circular_orbit_velocity(orbital_r) -> float
  - Takes in the orbital radius of an object
  - returns the speed at which it is traveling around the earth

"""
G = 6.6743 * math.pow(10,-11) # m3/(kg * s^2)
MASS_EARTH = 5.972 * math.pow(10,24) # kg
RADIUS_EARTH = 6371 # km

#returns the velocity in m/s for a object in circular orbit around the earth
def circular_orbit_velocity(orbital_r) -> float: 
  return math.sqrt((G * MASS_EARTH / orbital_r))

#returns the specific angular momentum as a vector (h = r X v)
def spec_ang_momentum_from_v_and_r(v, r) -> list:
  return vector_functions.cross_product(r,v)

#returns the wobble between two bodies in the given radius units
#If percent difference is true, will return a precent diff instead
#assumes m1 to be (0,0)
def wobble(m1, m2, r, percent_diff = False) -> float:
  wobble = (m2 * r)/(m2+m1) #compute center of mass
  #returns wobble or percent distance from m1 relative to r 
  return wobble if not percent_diff else wobble*100/r 

#returns approximately the orbital radius of an altitude
#alt should be in m from the earths surface
def altitude_to_orbital_radius_earth(alt, km=False) -> float:
  if km:
    return planetary_data.RADIUS_EARTH/1000 + alt 
  else:
    return planetary_data.RADIUS_EARTH + alt 

#returns a scalar state of distance and velocity. Distance measured from the central body
#returns None,None in the event that 1 + ecos(theta) is equal to 0 (div by 0 err)
#values are returned in m, m/s
def keplers_scalar_state(h, mu, e, theta) -> tuple: 

  if e < 0 or mu < 0: None, None

  if (1+e*math.cos(theta)) == 0 : return None, None 
  
  #Keplers equation r = h^2/mu^2 * (1/1+ecos(theta))
  radius = (math.pow(h,2)/mu)*(1/(1+e*math.cos(theta)))

  #Derived velocity equation
  velocity = (mu/h)*math.sqrt((1+e*(e+2*math.cos(theta))))

  return radius, velocity #return the value as a tuple

#expects si units (rp in m)
#From lecture slide hyperbolic trajectories summary of equations
#which cites curtis equations
def hyperbolic_state(e, rp, mu):

  a = rp/(e-1) # Curtis 2.105a
  ra = -a*(e+1) # Curtis 2.105b
  delta = 2 * math.asin(1/e) * (180/math.pi) #degrees 
  b = a * math.sqrt((e**2 - 1)) # Curtis 2.107
  theta_inf = math.acos(-1/e) * (180/math.pi) # Curtis 2.97 
  h = math.sqrt(mu * a * (e**2 - 1))
  beta = 180 - theta_inf
  v_inf = math.sqrt(mu/a) # Curtis 2.112 
  v_esc = math.sqrt(2*mu/rp)
  v_p = math.sqrt((v_esc**2) + (v_inf**2))
  epsilon = mu/(2*a)

  state = {
    'e'               : e,
    'rp (km)'          : rp,
    'ra (km)'          : ra, 
    'a (km)'           : a,
    'b (km)'           : b,
    'h (km^2/s)'       : h,
    'theta inf (deg)' : theta_inf,
    'Beta (deg)'      : beta,
    'delta (deg)'     : delta,
    'V_inf (km/s)'     : v_inf,
    'V_esc (km/s)'     : v_esc,
    'V_p (km/s)'     : v_p,
    'epsilon (MJ/kg)'  : epsilon
  }

  return state

#Returns a dictionary with all useful values of a defined orbit
#expects si units
#if km=True, pass in km values for a and mu
#else km=False, pass in si units
def orbital_state(a,e,mu, km=False) -> dict:
  r_p = a * (1-e) #rp = a(1-e)
  r_a = a * (1+e) #ra = a(1+e)
  if not km:
    T = (2*np.pi/np.sqrt(mu))*(a**(3/2)) # T = 2(pi) * sqrt(a^3/mu)
  else:
    #need Si units
    T = (2*np.pi/np.sqrt(mu*(10**9)))*((a*1000)**(3/2)) # T = 2(pi) * sqrt(a^3/mu)
  spec_e = -mu/(2*a) #-mu/(2*a)
  h = np.sqrt(a*mu*(1-np.pow(e,2))) # h = sqrt(mu*a*(1-e^2))
  b = a * np.sqrt(1-(np.pow(e,2)))
  v_p = h/r_p # vp = h/rp
  v_a = h/r_a # va = h/va

  if not km:
    state = {
      'e'             : e,
      'a (m)'         : a, 
      'r_p (m)'       : r_p, 
      'r_a (m)'       : r_a, 
      'v_p (m/s)'     : v_p,
      'v_a (m/s)'     : v_a,
      'b (m)'         : b,
      'h (m^2/s)'     : h,
      'T (s)'         : T,
      'spec_e (J/kg)' : spec_e
    }
  else:  
    state = {
      'e'              : e,
      'a (km)'         : a, 
      'r_p (km)'       : r_p, 
      'r_a (km)'       : r_a, 
      'v_p (km/s)'     : v_p,
      'v_a (km/s)'     : v_a,
      'b (km)'         : b,
      'h (km^2/s)'     : h,
      'T (s)'          : T,
      'spec_e (MJ/kg)' : spec_e,
      'Mu (km^3/s^2)'  : mu
    }
    
  return state

def print_state(orbital_state, label=None):
  if label:
    print(f'\n============ State of {label} ============')
  else:
    print('\n============ State ============')

  for state in orbital_state:
    print(f'{state} : {orbital_state[state]:.5g}')

def orbit_plot(h, mu, e, ax=None, color='b', units='km', label=None):
  if ax is None: #create plot if no plot parameter was given
        fig, ax = plt.subplots()
        ax.set_title('Orbital Plot')
        ax.set_ylabel(f'Distance ({units})')
        ax.set_xlabel(f'Distance ({units})')
        ax.plot(0, 0, 'r*', label='Planet CoM') #set focus to be 0,0
        

  theta = np.linspace(0, 2*np.pi,100)
  r = ((h**2)/mu)*(1/(1+(e*np.cos(theta)))) #radius plot

  #polar coords 
  x = r*np.cos(theta)
  y = r*np.sin(theta)

  ax.plot(x, y, color=color, label=label if label else f'e={e:.2f}')
  ax.set_aspect('equal', adjustable='datalim')
  formatter = ScalarFormatter(useMathText=True)
  formatter.set_scientific(True)
  formatter.set_powerlimits((-1, 1))

  ax.xaxis.set_major_formatter(formatter)
  ax.yaxis.set_major_formatter(formatter)
  ax.grid(True)

  ax.legend(*ax.get_legend_handles_labels(), loc="upper left")

  return ax #pass ax back to main 

#expects theta in degrees
#returns t in seconds
def t_from_theta(theta_deg, T, e):
  theta = theta_deg * math.pi / 180
  E = 2*math.atan((math.sqrt((1-e)/(1+e)))*math.tan(theta/2)) # Curtis 3.13b
  Me = E - (e * math.sin(E)) # Curtis 3.14
  t = (Me/(2*math.pi)) * T
  return t

def perifocal_plot(a, e, theta_deg, ax=None, color='b'):
    if ax is None: #create plot if no plot parameter was given
        fig, ax = plt.subplots()
  
    theta_rad = np.deg2rad(theta_deg)
    p_vec = np.array([np.cos(theta_rad), np.sin(theta_rad)]) #unit vector p
    q_vec = np.array([np.cos(theta_rad + np.pi/2), np.sin(theta_rad + np.pi/2)]) #unit vector q

    #create data points
    theta = np.linspace(0, 2*np.pi, 100)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * (p_vec[0]*np.cos(theta) + q_vec[0]*np.sin(theta))
    y = r * (p_vec[1]*np.cos(theta) + q_vec[1]*np.sin(theta))

    #Plot
    ax.plot(x, y, color=color)
    ax.plot(0, 0, 'r*') #set focus to be 0,0
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(True)

    return ax

#expects vectors in perifocal frame
def perifocal_delta_v(current_orbit, desired_orbit):
  vq = desired_orbit[1] - current_orbit[1]
  vp = desired_orbit[0] - current_orbit[0]
  return math.sqrt((vq**2) + (vp**2))

#returns in degrees
def flight_path_angle(v_t, v_r):
  return math.atan(v_r/v_t) * (180/math.pi)

#expects theta in rads
def compute_v_radial(h,mu,theta,e):
  return (mu/h)*(e*math.sin(theta))

#expects theta in rads
def compute_v_perpendicular(h,mu,theta,e):
  return (mu/h)*(1+e*math.cos(theta))

#Newton-Rhapson method Alg 3.1 
#returns in degrees
def theta_from_t(t, e, T):
  Me = t * (2*math.pi/T) # (3.8)
  E = Me + (e/2) if Me < math.pi else Me - (e/2)
  E_ratio = (E - e*math.sin(E)-Me)/(1 - e*math.cos(E))
  while(E_ratio > pow(10, -8)):
    E -= E_ratio 
    E_ratio = (E - e*math.sin(E)-Me)/(1 - e*math.cos(E));   

  e_ratio = math.sqrt((1+e)/(1-e))
  theta = (180/math.pi)*2*math.atan(e_ratio*math.tan(E/2))
  if theta < 0:
    return theta + 360
  else:
    return theta
  
#requires orbits to coaxial
#returns delta V
#expects units of km for length
def hohmann_transfer(init_state, final_state, mu, start_peri=True):

  if start_peri:
    #complete transfer at final apo 
    rp_xfer = init_state['r_p (km)']
    ra_xfer = final_state['r_a (km)']
  else:
    #complete transfer at final peri
    rp_xfer = final_state['r_p (km)']
    ra_xfer = init_state['r_a (km)']

  #find transfer orbit parameters
  a_xfer = (rp_xfer + ra_xfer) / 2
  e_xfer = (ra_xfer - rp_xfer) / (ra_xfer + rp_xfer)

  #get state velocities
  transfer_state = orbital_state(a_xfer, e_xfer, mu, km=True)
  va_xfer = transfer_state['v_a (km/s)']
  vp_xfer = transfer_state['v_p (km/s)']

  if start_peri:
    v_start_orbit = init_state['v_p (km/s)']  #burn at periapsis
    v_end_orbit   = final_state['v_a (km/s)']  #burn at apoapsis
    v_start_xfer  = vp_xfer#start transfer at peri
    v_end_xfer    = va_xfer
  else:
    v_start_orbit = init_state['v_a (km/s)']  #burn at apoapsis
    v_end_orbit   = final_state['v_p (km/s)'] #burn at periapsis
    v_start_xfer  = va_xfer#start transfer at apo
    v_end_xfer    = vp_xfer

  delta_v1 = abs(v_start_xfer - v_start_orbit) #burn 1
  delta_v2 = abs(v_end_orbit - v_end_xfer)     #burn 2 

  return (delta_v1, delta_v2, transfer_state) #full transfer magnitude

#Phasing manuever dtheta (expects signed angle)
#only works for circular orbits
#expects r & mu in si units
#more time = more fuel efficient
def phasing_maneuver(delta_theta_deg, r, mu, given_time, label = None):
  v_circ = math.sqrt(mu/r) #circular velo
  T = 2 * math.pi * math.sqrt((r**3)/mu) #period of a circular orbit
  n = (given_time/T)
  theta = delta_theta_deg * math.pi/180 # convert to rads
  T_xfer = T - (T*(theta)/(2* math.pi * n)) 
  a_xfer = (((T_xfer / (2 * math.pi)) ** 2) * mu) ** (1/3)

  #if positive start at peri (trailing)
  if a_xfer > r :
    e_xfer = 1  - (r/a_xfer)
  else:
    e_xfer = (r/a_xfer) - 1
  
  xfer_orbit = orbital_state(a_xfer, e_xfer, mu)

  if a_xfer < r:
    v_xfer = xfer_orbit['v_a (m/s)']
  else:
    v_xfer = xfer_orbit['v_p (m/s)']
  
  dv_total = 2 * abs(v_xfer-v_circ)

  if label:
    print(f'\n============= Transfer Orbit {label} =============')
  else:
    print('\n============= Transfer Orbit =============')

  print(f'Period of the Transfer orbit (s) {T_xfer:.1f}')
  print(f'Semi Major Axis of Transfer orbit (m) {a_xfer:.1f}')
  print(f'Eccentricity of Transfer orbit {e_xfer:.5f}')
  print(f'Total Needed Delta V in m/s: {dv_total:.4f}')

  return (a_xfer, e_xfer, dv_total)

#Does not work for times spent around periapsis because of the wrapping of angles 
#Finds the time spent between in an angular cone around apoapsis, expects range to be
#a singular number representing +/- a given angle from 180 
def time_spent_in_cone(range, e, T):
  t_to_range = t_from_theta(180-range,T,e)
  return (T/2 - t_to_range) * 2 

#need the parameter of the hyperbolic excess speed to solve
#assumes going to mars from earth 
#assumes km
def arrival_delta_v_mars(parking_alt):

  r_earth_sun = planetary_data.SUN_EARTH_KM
  r_mars_sun = planetary_data.SUN_MARS_KM
  mu_sun = planetary_data.MU_SUN_KM
  
  mu = planetary_data.MU_MARS_KM
  r_mars= planetary_data.RADIUS_MARS_KM + parking_alt #Find parking radius
  v_parking = np.sqrt(mu/r_mars) # Curtis 2.63
  v_esc = np.sqrt((2*mu)/r_mars) # Curtis 2.91

  # Curtis 8.35
  v_ab_a = np.sqrt((2*mu_sun*(r_earth_sun/r_mars_sun))/(r_earth_sun + r_mars_sun))
  v_inf = v_ab_a - planetary_data.V_MARS_SUN_KM # Curtis 8.51

  return abs(v_parking - np.sqrt((v_esc**2)+(v_inf**2))) # Curtis 8.40

#Total delta V needed to return from orbit of the parking alt
#Techically still orbiting just a 0 alt orbit
#assumes km
def return_from_parking_orbit_mars(parking_alt):
  
  mu = planetary_data.MU_MARS_KM
  r_parking = parking_alt + planetary_data.RADIUS_MARS_KM
  v_parking = np.sqrt(mu/r_parking) # 2.63 

  #find the states of these two orbits
  parking_state = orbital_state(r_parking, 0, mu, km=True)
  zero_orbit_state = orbital_state(planetary_data.RADIUS_MARS_KM, 0, mu, km=True)

  deltaV, transfer_orbit = hohmann_transfer(parking_state, zero_orbit_state, mu)

  return deltaV

def three_delta_V(parking_alt):
  deltaV1 = arrival_delta_v_mars(parking_alt)
  deltaV2_V3 = return_from_parking_orbit_mars(parking_alt)
  return deltaV1 + deltaV2_V3

def fourth_delta_v(parking_alt):
  #find the speed at which mars rotates
  angular_velo_mars = 2 * math.pi/planetary_data.SIDEREAL_ROTATION_MARS
  v_surface_mars = angular_velo_mars * planetary_data.RADIUS_MARS_KM
  
  #find the summed delta V
  mu = planetary_data.MU_MARS_KM
  delta_v_sum = three_delta_V(parking_alt)
  v_0_alt = np.sqrt(mu/planetary_data.RADIUS_MARS_KM) #parking speed

  return abs(v_0_alt - v_surface_mars) + delta_v_sum #total delta V 

#where r_inner and r_outer are relative to the sun
def get_v_ab_a(r_inner, r_outer):
  mu_sun = planetary_data.MU_SUN_KM
  # Curtis 8.35
  v_ab_a = np.sqrt((2*mu_sun*(r_inner/r_outer))/(r_inner + r_outer))
  return v_ab_a

def aiming_radius_from_parking_radius_mars(parking_alt):
  r_mars = planetary_data.RADIUS_MARS_KM
  mu = planetary_data.MU_MARS_KM
  parking_radius = parking_alt + r_mars

  v_ab_a = get_v_ab_a(
    planetary_data.SUN_EARTH_KM, planetary_data.SUN_MARS_KM)
  
  v_esc = np.sqrt((2*mu)/parking_radius) # Curtis 2.91
  v_inf = v_ab_a - planetary_data.V_MARS_SUN_KM # Curtis 8.51

  # at periapsis h = v_rp*r_rp
  h_hyper = parking_radius * np.sqrt((v_esc**2) + (v_inf**2))

  e = (h_hyper**2)/(mu * parking_radius) - 1 # Keplers 2nd 

  #get state of entry orbit
  hyper_state = hyperbolic_state(e, parking_radius, mu)

  return hyper_state['b (km)'] #return aiming radius

#v_inf
def get_v_inf_arrival_outer(r_inner, r_outer, planet_speed):
  v_ab_a = get_v_ab_a(r_inner, r_outer)
  return v_ab_a - planet_speed

#takes in ra and rp and returns the a & e value for the orbit
def a_e_from_rp_ra(rp,ra):
  e = (ra - rp)/(rp+ra) # Curtis 2.84
  a = (rp + ra)/2 # Definition of an orbit
  return (a,e)

#n is the mean angular velocity (2*pi)/T 
#returns the phase angle difference in degrees
#inner starts transfer at theta = 0 deg 
def hohmann_inner_outer_phase_difference(n, t):
  # Curtis 8.12
  return (math.pi - n*t) * 180/math.pi

def waiting_time_after_inner_outer(n_outer, n_inner, delta_phi):
  # Uses Curtis 8.15 
  return ((2*math.pi) - (delta_phi * math.pi/180))/(n_inner - n_outer)

def compute_delta_phi(n_inner, T_xfer):
  # Curtis 8.17
  return abs(2*(hohmann_inner_outer_phase_difference(n_inner, T_xfer/2))) 

def geocentric_3d_plot(omega, i, w, theta, a, e, deg=True):

  if deg:
      omega = np.deg2rad(omega)
      i = np.deg2rad(i)
      w = np.deg2rad(w)
      theta = np.deg2rad(theta)

  # Rotation matrices
  #Curtis 4.34
  R3_W = np.array([[ np.cos(omega),  np.sin(omega), 0],
                    [-np.sin(omega),  np.cos(omega), 0],
                    [0, 0, 1]])
  #Curtis 4.32
  R1_i = np.array([[1, 0, 0],
                    [0, np.cos(i), np.sin(i)],
                    [0, -np.sin(i), np.cos(i)]])
  #Curtis 4.34
  R3_w = np.array([[ np.cos(w),  np.sin(w), 0],
                    [-np.sin(w),  np.cos(w), 0],
                    [0, 0, 1]])

  #Rotation matrix - Curtis 4.49 
  Q_peri_geo = (R3_w @ R1_i @ R3_W).T 
  # the @ is syntactical sugar for matrix multiply! pretty awesome I just found
  # out about that while working on this 

  # Current radius vector
  r_mag = a*(1-e**2)/(1 + e*np.cos(theta))
  rc = r_mag * np.array([np.cos(theta), np.sin(theta), 0])
  rc_geo = Q_peri_geo @ rc 

  theta_deg = np.arange(0, 360)
  p = a*(1-e**2)/(1 + e*np.cos(np.deg2rad(theta_deg))) * np.cos(np.deg2rad(theta_deg))
  q = a*(1-e**2)/(1 + e*np.cos(np.deg2rad(theta_deg))) * np.sin(np.deg2rad(theta_deg))
  w_vec = np.zeros_like(theta_deg)
  r_peri = np.vstack((p, q, w_vec))
  r_geo = Q_peri_geo @ r_peri

  h = Q_peri_geo @ np.array([0,0,1])
  h = (a/2) * h # h_vec
  N = np.cross([0,0,1], h) # Node line
  rp = a*(1-e) * (Q_peri_geo @ np.array([1,0,0])) #rp 

  #equitorial planes
  x = RADIUS_EARTH * np.cos(np.deg2rad(theta_deg))
  y = RADIUS_EARTH * np.sin(np.deg2rad(theta_deg))
  z = np.zeros_like(theta_deg)
  eq_vec = np.vstack((x,y,z))

  # Plot
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.plot(r_geo[0], r_geo[1], r_geo[2], linewidth=1, label='Orbit')
  ax.plot(eq_vec[0], eq_vec[1], eq_vec[2],linewidth=1, label='Equitorial Plane', color = 'green')
  ax.quiver(0,0,0, rc_geo[0], rc_geo[1], rc_geo[2], color='red', linewidth=1, label='Current radius')
  ax.quiver(0,0,0, N[0], N[1], N[2], color='green', linestyle='-', label='Node line')
  ax.quiver(0,0,0, rp[0], rp[1], rp[2], color='orange', linestyle='-', label='Periapsis')
  ax.scatter([0],[0],[0], color='black', s=100, label='Earth')

  ax.set_xlabel('X [km]')
  ax.set_ylabel('Y [km]')
  ax.set_zlabel('Z [km]')
  ax.set_box_aspect([1,1,1])
  max_val = a*(1+e)
  ax.set_xlim([-max_val, max_val])
  ax.set_ylim([-max_val, max_val])
  ax.set_zlim([-max_val, max_val])
  ax.legend()
  plt.show()

#expects r and v vectors in the for i j k
#uses curtis alg 4.2 
def get_geocentric_orbital_params(r, v):

  #get magnitudes
  r_mag = vector_functions.magnitude(r)
  v_mag = vector_functions.magnitude(v)

  #Dot v with the unit vector of r to get velocity in radial dir
  v_r = (1/r_mag) * vector_functions.dot_product(r,v)
  
  #specific angular momentum vector 
  h = vector_functions.cross_product(r,v)
  h_mag = vector_functions.magnitude(h)

  # get inclination 
  i = math.acos(h[2]/h_mag) # 4.7

  #get node line vector k x h
  N = vector_functions.cross_product([0,0,1], h)
  N_mag = vector_functions.magnitude(N)

  #calculate right ascension 
  if (N[1] >= 0):
    omega = math.acos(N[0]/N_mag)
  else:
    omega = math.pi * 2 - math.acos(N[0]/N_mag)

  mu = planetary_data.MU_EARTH_KM

  #get e vector
  vec1 = [((v_mag**2)- (mu/r_mag)) * component for component in r]
  vec2 = [(r_mag * v_r) * component for component in v]
  vec3 = vector_functions.subtraction(vec1, vec2)
  e = [(1/mu) * component for component in vec3]

  #get eccentricity
  e_mag = vector_functions.magnitude(e)

  #arument of perigee
  if(e[2] >= 0):
    w = math.acos(vector_functions.dot_product(N,e)/(N_mag * e_mag))
  else:
    w = math.pi * 2 - math.acos(vector_functions.dot_product(N,e)/(N_mag * e_mag))
  
  #get the true anomaly 
  e_unit = vector_functions.unit_vector(e)
  r_unit = vector_functions.unit_vector(r)

  if(v_r >= 0):
    theta = math.acos(vector_functions.dot_product(e_unit, r_unit))
  else:
    theta = math.pi * 2 - math.acos(vector_functions.dot_product(e_unit, r_unit))

  state = {
    'h (km^2/s)'  : h_mag, 
    'omega (deg)' : np.rad2deg(omega),
    'i (deg)'     : np.rad2deg(i),
    'w (deg)'     : np.rad2deg(w),
    'e'           : e_mag,
    'theta (deg)' : np.rad2deg(theta)
  }

  print(
    f'Node Vector {fmt_list(N)} Km^2/s \n'
    f'H vector {fmt_list(h)} km^2/s\n'
    f'e vector {fmt_list(e)}'
  )

  return state

def fmt_list(lst, precision=".5g"):
  """Return a string containing the list with each element formatted."""
  return "[" + ", ".join(f"{x:{precision}}" for x in lst) + "]"

def alg_5_4(A, Ar, a, ar, p, pr, lat, sidr, H, deg=True):

  if deg:
    A  = vector_functions.deg2rad(A)
    Ar = vector_functions.deg2rad(Ar)
    a  = vector_functions.deg2rad(a)
    ar = vector_functions.deg2rad(ar)
    lat = vector_functions.deg2rad(lat)
    sidr = vector_functions.deg2rad(sidr)

  f = planetary_data.EARTH_FLATTENING_FACTOR
  Re = planetary_data.RADIUS_EARTH_KM

  big_term = (Re/math.sqrt(1-(2*f - (f**2))*(math.sin(lat)**2)))

  # Curtis 5.56
  R = [
    (big_term + H)*math.cos(lat)*math.cos(sidr),
    (big_term + H)*math.cos(lat)*math.sin(sidr),
    ((big_term*(1-f)**2) + H)*math.sin(lat)
  ]

  #5.83a
  delta = math.asin(math.cos(lat)*math.cos(A)*math.cos(a) + math.sin(lat)*math.sin(a))

  #5.83b
  temp = math.acos((math.cos(lat)*math.sin(a) - math.sin(lat)*math.cos(A)*math.cos(a))/math.cos(delta))
  h = math.pi * 2 - temp if A < math.pi else temp

  #5.83c
  alpha = sidr - h 

  #5.68 & 5.71
  p_unit = [
    math.cos(delta)*math.cos(alpha),
    math.cos(delta)*math.sin(alpha),
    math.sin(delta)
  ]

  p_vec = [ele * p for ele in p_unit]

  #5.63
  r = vector_functions.addition(R, p_vec)

  #2.67
  omega = [0,0,planetary_data.ANG_VEL_EARTH]
  R_dot = vector_functions.cross_product(omega, R)

  #5.84
  delta_dot = (1/math.cos(delta))*(
    (-Ar*math.cos(lat)*math.sin(A)*math.cos(a)) + 
    ar * (math.sin(lat)*math.cos(a) - math.cos(lat)*math.cos(A)*math.sin(a)))
  
  #5.85
  alpha_dot = planetary_data.ANG_VEL_EARTH + (
    (Ar*math.cos(A)*math.cos(a) - ar*math.sin(A)*math.sin(a) + delta_dot*math.sin(A)*math.cos(a)*math.tan(delta))/
    (math.cos(lat)*math.sin(a) - math.sin(lat)* math.cos(A)* math.cos(a))
  )

  # 5.69 
  pr_unit= [
    (-alpha_dot*math.sin(alpha)*math.cos(delta) - delta_dot * math.cos(alpha) * math.sin(delta)),
    (alpha_dot*math.cos(alpha)*math.cos(delta) - delta_dot * math.sin(alpha) * math.sin(delta)),
    (delta_dot * math.cos(delta))
  ]

  p1 = [pr * ele for ele in p_unit]
  p2 = [p * ele for ele in pr_unit]

  v = vector_functions.addition(R_dot, vector_functions.addition(p1, p2))

  return (r, v)

#expects si units
def arrival_depart_surface_planetary(T_sidereal, alt_parking, mu, r_planet, departure):

  #Delta V to get to 0 alt orbit
  Vsurf = ((2*math.pi)/T_sidereal) * r_planet #surface velo w*r
  Vzero = math.sqrt(mu/r_planet) # Curtis 2.63
  dVsurf = abs(Vzero - Vsurf)

  #Get states
  r_parking = alt_parking + r_planet

  if departure:
    initial_state = orbital_state(r_planet, 0, mu, km=True)
    final_state = orbital_state(r_parking, 0, mu, km=True)
  else:
    initial_state = orbital_state(r_parking, 0, mu, km=True)
    final_state = orbital_state(r_planet, 0, mu, km=True)

  #Hohmann Transfer
  if departure:
    dv1, dv2, transfer_orbit = hohmann_transfer(
      initial_state, final_state, mu, start_peri=True
    )

    #       Departure, Arrival
    return (dv1 + dVsurf, dv2)
  else:
    dv1, dv2, transfer_orbit = hohmann_transfer(
      initial_state, final_state, mu, start_peri=False
    )
    #       Departure, Arrival
    return (dv1, dv2 + dVsurf)
  
def arrival_depart_interplanetary(
    Rarr, Rdept, Vdept, Varr, Vpark_dept, Vpark_arr, Vesc_dept, Vesc_arr):
  mu = planetary_data.MU_SUN_KM

  Vab_dept = math.sqrt(2*mu*(Rarr/Rdept)/(Rarr + Rdept))
  Vinf_dept = Vab_dept - Vdept
  Vhyper_peri_dept = math.sqrt((Vesc_dept**2) + (Vinf_dept**2)) 

  Vab_arr = math.sqrt(2*mu*(Rdept/Rarr)/(Rarr + Rdept))
  Vinf_arr = Vab_arr - Varr
  Vhyper_peri_arr = math.sqrt((Vesc_arr**2) + (Vinf_arr**2)) 
  
  #           DVDepatrure                        DVArrival              
  return abs(Vhyper_peri_dept - Vpark_dept), abs(Vhyper_peri_arr - Vpark_arr)

def get_v_esc(mu, r_peri):
  return math.sqrt(2*mu/r_peri) # Curtis 2.91

def get_circular_orbit_speed(mu, r):
  return math.sqrt(mu/r) # Curtis 2.63

def get_mass_prop_needed_for_dv(dv, c, m_init):
  m_final = m_init * pow(math.e, -(dv/c)) # Rocket Equation

  #          propellant mass used     Final Mass
  return (    m_init - m_final     ,    m_final    )

def get_needed_prop_mass(m_final, dv, c):
  return m_final * (pow(math.e, (dv/c)) - 1) # Rocket Equation

# Assumes a plane change maintaining the speed before and after 
# Uses Curtis 6.28
def inclination_change(v, inclination, deg=True):

  if deg:
    inclination = np.deg2rad(inclination)
  
  return 2 * v * math.sin(inclination/2)

# Get the phase difference L11, Interplanetary Transfers 
# returns in radians
def get_phase_difference_in_out(Txfer, n):
  t = Txfer/2
  phi_dept = math.pi - (n*t)
  return phi_dept

# Get the period
def get_period(a,mu):
  T = (2*np.pi/np.sqrt(mu))*(a**(3/2)) # Curtis 2.83
  return T

def get_eccentricity_turn_angle(delta, deg=True):

  if deg:
    delta = np.deg2rad(delta)

  # Curtis 2.100
  return 1/(math.sin(delta/2))

  



