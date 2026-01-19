""" Q 2 """

import math

print('P2')

R = 50 
f = 10 * 10**3
w = math.pi * 2 * f

v_out_over_v_in = (1/5)

L = (R/w) * math.sqrt((1/(v_out_over_v_in**2)) - 1)

print(f' Inductance {L * (10**3)} mH')

""" Q 4 """

print("\nProblem 4")

R = 200 
L = 0.400 
C = 6 * pow(10,-6)

w = 250
Vp = 30 

Xc = 1/(w*C)
Xl = w*L
Z = math.sqrt((R**2) + (Xl - Xc)**2) #impedance

print(f'Impedance {Z} Ohms') 

Ip = Vp/Z #A 

print(f'Current Amp {Ip * (10**3)} mA')

phi = -(180/math.pi) * math.acos(R/Z) #deg

if(Xl > Xc): print("Voltage Leads")
else: print('Voltage Lags')

print(f'Phase angle {phi} deg')

w_res = math.sqrt(1/(L*C)) 
f_res = w_res/(2*math.pi)

print(f'resonant freq {f_res} Hz')

I_res = Vp/R

print(f'Curr at res {I_res * (10**3)} mA')

""" It's very unintuitive as to why you can measure voltages across the capacitor and inductor to be greater than the source voltage. The reason is that these voltage phasors (inductor and capacitor) are out of phase with each other, to be precise, 180 degrees out of phase. This results in these voltages cancelling each other out, and so at the end of the day, Kirchhoff's loop rule still holds. This is also why at resonance we see the peak current, at this point the magnitudes of the phasors are equal, and so they fully cancel each other out, and we see the voltage of the source go across the resistor (the maximum current). """

""" Q 11 """

print("\nQ11")
Vs = 15 #secondary
Vp = 120 #primary

winds_ratio = (Vs/Vp)

print(f'Winds ratio {winds_ratio}')

Vpeak = Vs * math.sqrt(2)

print(f'Peak {Vpeak} V')

Vpeak_peak = Vpeak * 2 

print(f'Peak to Peak {Vpeak_peak} V')

R = 30

P_rms = (Vs**2)/R

print(f'Power rms Dissipated {P_rms} W')

P_peak = (Vpeak**2)/R

print(f'Power peak dissipated {P_peak} W')

Is_peak = Vpeak/R

Ip_peak = Is_peak * (winds_ratio)

print(f'Peak current in the VAC source {Ip_peak * (10**3)} mA')


