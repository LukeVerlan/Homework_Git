from matplotlib import pyplot as plt
import numpy as np

inclination = np.array([
  10, 20, 30, 40, 50, 60, 70
])

Vout = np.array([
  0.41, 0.80, 1.19, 1.56, 1.78, 2.05, 2.26
])

fig, ax = plt.subplots()

ax.plot(np.sin(np.deg2rad(inclination)), Vout, color='g')
ax.set_title('V Vs Sin(θ)')
ax.set_xlabel('Sin(θ)')
ax.set_ylabel('V')

k = Vout / np.sin(np.deg2rad(inclination))

kavg = sum(k)/len(k)

print(f' k_avg: {kavg:.5g}')

T0 = 40
w0 = 1000
Vs = 25

w_range = np.arange(0, 1000)
T = T0 * (1 - (w_range/w0))
T_l = 0.5*(w_range*2*np.pi)/60 #convert to rad/s

fig, ax1 = plt.subplots()
ax1.plot(w_range, T, color='b', label='Unloaded')
ax1.plot(w_range, T_l, color='r', label='loaded')

ax1.set_title('T vs RPM')
ax1.set_xlabel('RPM')
ax1.set_ylabel('T N*m')

# plt.show()

import math

omega = 1/((0.5/T0)+(1/w0))
T = omega * 0.5 
P = T * (math.pi * 2 * omega/60)
I = P/Vs

print(f'omega in rpm {omega:.5g}')
print(f'I in amps {I:.5g}')

Vs_ratio = 40/Vs
w0 *= Vs_ratio
T0 *= Vs_ratio
omega = 1/((0.5/T0)+(1/w0))
T = omega * 0.5
print(f'omega in rpm {omega:.5g}')
P = T * (omega * 2 *math.pi/60)
I = P/40
print(f'I in amps {I:.5g}')
