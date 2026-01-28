import elasticity_tools as et 


def main(): 
    problem_1_3()

def problem_1_3():
    states = [
        (54, 30, 5),
        (30, 54, -5),
        (-60, -36, 5),
        (30, -50, 30)
    ]

    for i in range(len(states)):
        state = states[i]
        stress_state = et.plane_stress_state(state[0],state[1],state[2],stress_units="N/mm^2")
        et.print_state(stress_state, label = f"State No. {i+1}")
        et.mohrs_circle_plot(stress_state, title=f"State No. {i+1}")

def problem_1_1():
    sig_x = 80
    sig_y = 0
    tau_xy = 45

    state = et.plane_stress_state(
        sig_x,sig_y, tau_xy, stress_units="N/mm^2"
    )

    et.print_state(state)

if __name__ == "__main__":
    main()