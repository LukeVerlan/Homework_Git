import math
import numpy
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

    OFFSET = stress_state['tau_max'][0]/10
    
    if stress_state['sig_y'][0] > stress_state['sig_x'][0]:
        sx_offset = -OFFSET
        sy_offset = OFFSET
    else: 
        sx_offset = OFFSET
        sy_offset = -OFFSET

    t_offset = OFFSET

    fig, ax = plt.subplots(figsize=(6, 6)) 

    stress_units = stress_state['sig_1'][1]

    ax.set_title(title, fontweight='bold')
    ax.set_ylabel(f"Shear Stress ({sym.tau}), {stress_units}")
    ax.set_xlabel(f"Normal Stress ({sym.sigma}), {stress_units}")
    
    ax.invert_yaxis()

    ax.axhline(0, color='black', linewidth=1.5, zorder=2)
    ax.axvline(0, color='black', linewidth=1.5, zorder=2)
    ax.grid(True, linestyle='--', alpha=0.6, zorder=1)

    # Circle Base
    avg_stress = stress_state['sig_avg'][0]
    radius = stress_state['tau_max'][0]
    
    circle = Circle((avg_stress, 0), radius, color='blue', fill=False, linewidth=2.5, zorder=3)
    ax.add_patch(circle)

    sig_coords = [ stress_state['sig_y'][0], stress_state['sig_x'][0]]
    tau_coords = [-stress_state['tau_xy'][0], stress_state['tau_xy'][0]]
    ax.plot(sig_coords, tau_coords, color='red', linestyle='-', linewidth=2, zorder=4)

    padding = radius * 0.2
    ax.set_xlim(avg_stress - radius - padding, avg_stress + radius + padding)
    ax.set_ylim(radius + padding, -radius - padding)

    # Key points
    ax.plot(stress_state['sig_x'][0], stress_state['tau_xy'][0], marker='o', color='red')
    ax.annotate(f'{sym.sigma}x', xy=((stress_state['sig_x'][0], stress_state['tau_xy'][0])), 
                xytext=((stress_state['sig_x'][0]+ sx_offset, stress_state['tau_xy'][0]+ t_offset)))

    ax.plot( stress_state['sig_y'][0], -stress_state['tau_xy'][0], marker='o', color='red')
    ax.annotate(f'{sym.sigma}y', xy=((stress_state['sig_y'][0], -stress_state['tau_xy'][0])), 
                xytext=((stress_state['sig_y'][0]+ sy_offset, -stress_state['tau_xy'][0]- t_offset)))
    
    ax.plot( stress_state['sig_avg'][0], 0, marker='o', color='red')
    ax.annotate(f'{sym.sigma}avg', xy=((stress_state['sig_avg'][0], 0)), 
                xytext=((stress_state['sig_avg'][0], sy_offset)))
    
    ax.plot( stress_state['sig_avg'][0], stress_state['tau_max'][0], marker='o', color='red')
    ax.annotate(f'{sym.tau}max', xy=((stress_state['sig_avg'][0], stress_state['tau_max'][0])), 
                xytext=((stress_state['sig_avg'][0], stress_state['tau_max'][0] + t_offset)))

    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()

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

# Expects E in si units
def plane_strain_to_stress(strain_state, E, v):
    if strain_state['ex'][1] != sym.mu: fac = 1
    else:                               fac = pow(10, -6)
    G = E / (2*(1+v))
    
    # returns sig1
    def strain_to_stess(e1, e2): return E * fac * (e1 + v * e2) / (1 - (v**2))

    sigy = strain_to_stess(
        strain_state['ey'][0], strain_state['ex'][0]
    ) * pow(10, -6)

    sigx = strain_to_stess(
        strain_state['ex'][0], strain_state['ey'][0]
    ) * pow(10, -6)

    t_xy = strain_state['y_xy'][0] * fac * G * pow(10, -6) 

    return plane_stress_state(
        sigx, sigy, t_xy, 'MPa'
    )

def mohrs_strain_circle_plot(strain_state, title=""):

    OFFSET = strain_state['y_max'][0]/20
    
    if strain_state['ey'][0] > strain_state['ex'][0]:
        sx_offset = -OFFSET
        sy_offset = OFFSET
    else: 
        sx_offset = OFFSET
        sy_offset = -OFFSET

    t_offset = OFFSET

    fig, ax = plt.subplots(figsize=(6, 6)) 

    strain_units = strain_state['ex'][1]

    ax.set_title(title, fontweight='bold')
    ax.set_ylabel(f"Shear Strain ({sym.gamma}/2), {strain_units}")
    ax.set_xlabel(f"Normal Strain ({sym.epsilon}), {strain_units}")
    
    ax.invert_yaxis()

    ax.axhline(0, color='black', linewidth=1.5, zorder=2)
    ax.axvline(0, color='black', linewidth=1.5, zorder=2)
    ax.grid(True, linestyle='--', alpha=0.6, zorder=1)

    # Circle Base
    avg_strain = strain_state['eavg'][0]
    radius = strain_state['y_max'][0]/2
    
    circle = Circle((avg_strain, 0), radius, color='blue', fill=False, linewidth=2.5, zorder=3)
    ax.add_patch(circle)

    ey = strain_state['ey'][0]
    ex = strain_state['ex'][0]
    yxy = strain_state['y_xy'][0]/2

    sig_coords = [ ey, ex ]
    tau_coords = [-yxy, yxy]
    ax.plot(sig_coords, tau_coords, color='red', linestyle='-', linewidth=2, zorder=4)

    padding = radius * 0.2
    ax.set_xlim(avg_strain - radius - padding,  avg_strain + radius + padding)
    ax.set_ylim(radius + padding, -radius - padding)

    # Key points
    ax.plot(ex, yxy, marker='o', color='red')
    ax.annotate(f'{sym.epsilon}x', xy=((ex, yxy)), 
                xytext=((ex + sx_offset, yxy + t_offset)))

    ax.plot( ey , -yxy, marker='o', color='red')
    ax.annotate(f'{sym.epsilon}y', xy=((ey, -yxy)), 
                xytext=((ey + sy_offset, -yxy - t_offset)))
    
    ax.plot( avg_strain, 0, marker='o', color='red')
    ax.annotate(f'{sym.epsilon}avg', xy=((avg_strain, 0)), 
                xytext=((avg_strain, sy_offset)))
    
    ax.plot( avg_strain, radius, marker='o', color='red')
    ax.annotate(f'{sym.gamma}max/2', xy=((avg_strain, radius)), 
                xytext=((avg_strain, radius + t_offset)))

    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.show()