import math

m = 10
kp = 6
ks = 3/2
bp = 8
bs = 2
k  = 3
b  = 4

options = [
#   k, b, option letter
    (kp, b, 'a'),
    (ks, b, 'b'),
    (k, bp, 'c'),
    (k, bs, 'd'),
    (kp, bp, 'e'),
    (ks, bp, 'f'),
    (kp, bs, 'g'),
    (ks, bs, 'h')
]

def condition(K, B): return B/math.sqrt(2*K*m)

print("Successful Options:")
for option in options: 
    if condition(option[0], option[1]) >= 1: print(f"   {option[2]}")
