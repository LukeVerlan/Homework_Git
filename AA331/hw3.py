import elasticity_tools as et 
import math
import numpy as np
from matplotlib import pyplot as plt

import sys 
sys.path.append('../')

import symbols as sym


def problem_3_b():
    ex = 1100 * pow(10,-6)
    ey = -300 * pow(10,-6)
    e3 = -400 * pow(10,-6)

    E = 200 * (10**9)
    v = 0.3

    d = 75 * pow(10,-3)

    # Equation (1.31)
    y_xy = (e3 - ex*(math.cos(math.pi/4)**2) - ey * (math.sin(math.pi/4)**2)) / \
           (math.cos(math.pi/4) * math.sin(math.pi/4))
    
    strain_state = et.plane_strain_state(
       ex, ey, y_xy
    )

    # Strain state
    et.print_state(strain_state, label='Strain State')

    stress_state = et.plane_strain_to_stress(
        strain_state, E, v
    )

    et.print_state(stress_state, label='Stress State')

    et.mohrs_strain_circle_plot(
        strain_state, title='Mohrs Circle of Strain'
    )

    et.print_axial_torsional_loading(stress_state, d/2, pow(10,6))

    """ PART C """
    P = 980.7 * pow(10,3) #KN
    T = 10.2 * pow(10,3) #KN*m 

    angles = np.linspace(0, 90, 100)

    alpha = 45
    beta  = 45
    ex # axial direction, eII
    ey # longitudinal direction eI
    y_xy # the shear stress term

    # Equation (1.31)
    def rotation(angle):
        return (1/2)*(ey+ex) + (1/2)*(ey-ex) * np.cos(2*(np.deg2rad(angle))) + \
               (y_xy/2) * np.sin(2*(np.deg2rad(angle)))

    e2 = rotation(angles)
    e3 = rotation(angles + alpha)
    e1 = rotation(angles + alpha + beta)

    fig, ax = plt.subplots()
    
    ax.set_title('e1, e2, e3 Vs rotation angle')
    ax.set_xlabel('Degrees')
    ax.set_ylabel(f'Strain {sym.epsilon}')
    ax.plot(angles, e1, label='e1', color='r')
    ax.plot(angles, e2, label='e2', color='b')
    ax.plot(angles, e3, label='e3', color='g')
    ax.legend(loc='upper left')
    
    plt.show()

def problem_4():
    E = 73 * pow(10,9)
    v = 0.3

    ex = 3657 #mu
    ey = -1245 
    e3 = 956
    theta = 45
    
    yxy = et.get_shear_strain_from_rosette(
        e3, theta, ex, ey
    )

    strain_state = et.plane_strain_state(
        ex, ey, yxy, strain_units=sym.mu
    )

    et.print_state(
        strain_state, label='Strain State'
    )

    et.mohrs_strain_circle_plot(
        strain_state, title='Strain Circle'
    )

    stress_state = et.plane_strain_to_stress(
        strain_state, E, v
    )

    et.print_state(
        stress_state, label='Stress State'
    )

def main():
    problem_4()

if __name__ == '__main__':
    main()