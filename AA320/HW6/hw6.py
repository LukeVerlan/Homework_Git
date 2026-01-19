import math
from matplotlib import pyplot as plt 

""" Question 1 """

#Reff = Rf/G

Rf = 100 * pow(10,3)

G = 10**6

print(f'Reff = {Rf/G:.5g} ohms')

""" Question 2 """

Rf # same as last

Responstivity = 0.3 #uA/uW

Vout = 0.2 # V 

#adjust for uA
Pl = (Vout*(10**6))/(Rf*Responstivity)

print(f'Power out = {Pl:.5g} uA/uW')

""" Question 3 """

R = 20 * pow(10,3)
C = 100 * pow(10, -6)

tau = R*C

