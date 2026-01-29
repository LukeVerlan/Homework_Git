import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

import sys 
sys.path.append('../')

import symbols as sym


# =======================================
#               STRESS 
# =======================================

def plane_stress_state(sig_x, sig_y, tau_xy, stress_units="", angle_units="deg") -> dict:
    sig_avg = (sig_x + sig_y)/2
    radius  = (1/2) * math.sqrt((sig_x - sig_y)**2 + 4 * (tau_xy**2))
    sig_1   = sig_avg + radius
    sig_2   = sig_avg - radius
    theta_p = (1/2)*math.atan(2 * tau_xy/(sig_x - sig_y))
    if angle_units == "deg": theta_p *= (180/math.pi)

    state = {
        'sig_x'         : (sig_x, stress_units),
        'sig_y'         : (sig_y, stress_units),
        'tau_xy'        : (tau_xy, stress_units),
        'sig_avg'       : (sig_avg, stress_units),
        'tau_max'       : (radius, stress_units),
        'sig_1'         : (sig_1, stress_units),
        'sig_2'         : (sig_2, stress_units),
        'theta_p'       : (theta_p, angle_units),
        'theta_p_mohrs' : (2*theta_p, angle_units)
    }

    return state

def print_state(state, label=None):
    if label : print(f"================= {label} ==================")
    else     : print("============================================")

    for key in state: print(f"    {key}: {state[key][0]:.4g}, {state[key][1]}")

    print("============================================")

def mohrs_circle_plot(stress_state, title=""):

    circle_plot(
        stress_state['sig_x'][0], stress_state['sig_y'][0], stress_state['tau_xy'][0],
        stress_state['tau_max'][0], stress_state['sig_avg'][0], title=title,
        units=stress_state['sig_x'][1], xlabel=f'Normal stress ({sym.sigma})',
        ylabel=f'Shear Stress {sym.tau}', xsym=sym.sigma, ysym=sym.tau
    )

# Assumes a circular rod J = pi * r**4 / 2
# Assumes x is axial
# Stress Scale factor is for non si units, ie for MPa it should be 10^6
def print_axial_torsional_loading(stress_state, r, stress_scale_factor):
    
    P = stress_scale_factor * abs(stress_state['sig_x'][0]) * math.pi * (r**2)
    T = stress_scale_factor * abs(stress_state['tau_xy'][0]) * math.pi * ((r**4)/2) / r

    print(f'Axial Load in N {P:.5g} Torque in N*m {T:.5g}')


# =======================================
#               STRAIN
# =======================================

# Assumes 2 gauges are perpendicular, and the 3rd is an angle theta from the x axis
def get_shear_strain_from_rosette(e3, angle, ex, ey, deg=True):
    if deg: angle = np.deg2rad(angle)
    return 2 * (e3 - ex * (math.cos(angle)**2) - ey * (math.sin(angle)**2)) / (math.sin(2*angle))



def plane_strain_state(ex, ey, y_xy, strain_units='', angle_units='deg'):
    eavg = (ex + ey)/2
    R = (1/2)*math.sqrt((ex-ey)**2 + y_xy**2)
    e1 = eavg + R
    e2 = eavg - R
    theta_p = (1/2) * math.atan2(y_xy,ex - ey)
    y_max = 2*R 
    if angle_units == "deg": theta_p *= (180/math.pi)
    # else radians

    state = {
        'ex'            : (ex, strain_units),
        'ey'            : (ey, strain_units),
        'y_xy'          : (y_xy, strain_units),
        'eavg'          : (eavg, strain_units),
        'y_max'         : (y_max, strain_units),
        'e1'            : (e1, strain_units),
        'e2'            : (e2, strain_units),
        'theta_p'       : (theta_p, angle_units),
        'theta_p_mohrs' : (2*theta_p, angle_units)
    }

    return state

# Expects SI units
def plane_strain_to_stress(strain_state, E, v):

    if strain_state['ex'][1] == sym.mu:
        fac = 1e-6
    else:
        fac = 1.0
        
    G = E / (2 * (1 + v))
    
    def calc_normal_stress(e1, e2): return (E / (1 - v**2)) * (e1 * fac + v * (e2 * fac))

    sigx_pa = calc_normal_stress(strain_state['ex'][0], strain_state['ey'][0])
    sigy_pa = calc_normal_stress(strain_state['ey'][0], strain_state['ex'][0])

    txy_pa = G * (strain_state['y_xy'][0] * fac)

    return plane_stress_state(
        sigx_pa / 1e6, 
        sigy_pa / 1e6, 
        txy_pa / 1e6, 
        'MPa'
    )

def mohrs_strain_circle_plot(strain_state, title=""):

    circle_plot(
        strain_state['ex'][0], strain_state['ey'][0], strain_state['y_xy'][0]/2,
        strain_state['y_max'][0]/2, strain_state['eavg'][0], title=title,
        units=strain_state['ex'][1], xlabel=f'Normal Strain ({sym.epsilon})',
        ylabel=f'Shear Strain {sym.gamma}/2', xsym=sym.epsilon, ysym=sym.gamma
    )

def circle_plot(x, y, xy, r, avg, 
                title='', units='',
                xlabel='', ylabel='',
                xsym='', ysym=''):
    OFFSET = r/10
    
    if y > x:
        x_offset = -OFFSET
        y_offset = OFFSET
    else: 
        x_offset = OFFSET
        y_offset = -OFFSET

    t_offset = OFFSET

    fig, ax = plt.subplots(figsize=(6, 6)) 

    ax.set_title(title, fontweight='bold')
    ax.set_ylabel(f"{ylabel}, {units}")
    ax.set_xlabel(f'{xlabel}, {units}')
    
    ax.invert_yaxis()

    ax.axhline(0, color='black', linewidth=1.5, zorder=2)
    ax.axvline(0, color='black', linewidth=1.5, zorder=2)
    ax.grid(True, linestyle='--', alpha=0.6, zorder=1)
    
    circle = Circle((avg, 0), r, color='blue', fill=False, linewidth=2.5, zorder=3)
    ax.add_patch(circle)

    xcoords = [ y, x ]
    ycoords = [-xy, xy]
    ax.plot(xcoords, ycoords, color='red', linestyle='-', linewidth=2, zorder=4)

    padding = r * 0.2
    ax.set_xlim(avg - r - padding,  avg + r + padding)
    ax.set_ylim(r + padding, -r - padding)

    # Key points
    ax.plot(x, xy, marker='o', color='red')
    ax.annotate(f'{xsym}x', xy=((x, xy)), 
                xytext=((x + x_offset, xy + t_offset)))

    ax.plot( y , -xy, marker='o', color='red')
    ax.annotate(f'{xsym}y', xy=((y, -xy)), 
                xytext=((y + y_offset, -xy - t_offset)))
    
    ax.plot( avg, 0, marker='o', color='red')
    ax.annotate(f'{xsym}avg', xy=((avg, 0)), 
                xytext=((avg, y_offset)))
    
    ax.plot( avg, r, marker='o', color='red')
    ax.annotate(f'{ysym}max/2', xy=((avg, r)), 
                xytext=((avg, r + t_offset)))

    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()