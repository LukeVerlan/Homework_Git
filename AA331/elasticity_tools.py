import math
import numpy
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

import sys 
sys.path.append('../')

import symbols as sym

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
