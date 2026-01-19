"""
This file contains vector functions to be used in other files
Functions here:
  - dot_product(vector a, vector b) -> float 
		- Vectors must be lists of the same length and have same var indexes

  - cross_product(vector a, vector b) -> list 
		- Vectors must be lists of length 3 and have same var indexes
		- follows the format of (a X b)

	- scalar_triple_product(a,b,c) -> float
		- Vectors must be lists of the same length and have same var indexes
		- follows the format of (a * (b X c))

	- vector_triple_product(a,b,c) -> list
		- Vectors must be lists of length 3 and have same var indexes
		- follows the format of (a X (b X c))

"""

import math

def dot_product(vector_a, vector_b) -> float: #returns a float(scalar placeholder) 
	if len(vector_a) != len(vector_b): raise IndexError # Vectors must be same length
	product = 0 #init sum
	for var in range(len(vector_a)): product += vector_a[var] * vector_b[var] #add the product of all vars
	return product	#return the value

#indexes
#	- 0 : x 
# - 1 : y
# - 2 : z
def cross_product(vector_a, vector_b) -> list: #returns a list vector
	if len(vector_a) != len(vector_b) or len(vector_b) != 3 : raise IndexError #Vectors must be of R_3 and same length
	return 	[
							vector_a[1]*vector_b[2] - vector_a[2]*vector_b[1], #(a_y*b_z - a_z*b_y)
							vector_b[0]*vector_a[2] - vector_a[0]*vector_b[2], #(b_x*a_z - a_x*b_z)
							vector_a[0]*vector_b[1] - vector_b[0]*vector_a[1]	 #(a_x*b_y - b_x*a_y)
					]

#returns a scalar
#performs a dot product operation of A on the resulting cross product of B X C 
def scalar_triple_product(a,b,c) -> float: return dot_product(a, cross_product(b,c))

#returns a vector
#performs a cross product operation of A on the resulting cross product of B X C 
def vector_triple_product(a,b,c) -> list: return cross_product(a, cross_product(b,c))

#returns a scalar, takes in a vector (list) of any size
def magnitude(v) -> float:
	inside_sqrt = 0
	for var in v : inside_sqrt += math.pow(var, 2) #square each value
	return math.sqrt(inside_sqrt) #take sqrt & return

#subtracts v_subtractor from v, expects both vectors to be of the same length
def subtraction(v, v_subtractor):
	result = []
	for i in range(len(v)):
		result.append(v[i] - v_subtractor[i])
	return result

#adds two vectors
def addition(v1, v2):
  v_res = []
  for i in range(len(v1)): v_res.append(v1[i] + v2[i])
  return v_res

def matrix_transpose(matrix):
    return [list(row) for row in zip(*matrix)]

def matrix_multiply(A, B):
    # Sizes
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])

    # Check dimension compatibility
    if cols_A != rows_B:
        raise ValueError("Incompatible dimensions: A columns must equal B rows.")

    # Create output matrix filled with zeros
    C = [[0] * cols_B for _ in range(rows_A)]

    # Multiply
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):  # same as rows_B
                C[i][j] += A[i][k] * B[k][j]

    return C

#Converts degrees to radians
def deg2rad(deg):
  return deg*math.pi/180; 

#converts radians to degrees
def rad2deg(rad):
  return rad * 180/math.pi

# unit_vector
def unit_vector(v):
  mag = magnitude(v)
  return [(1/mag) * component for component in v]


def print_vector_formatted(l, vec_units, units='', label=None):
	if label: print(f'{label}: ', end=' ')
	for i in range(len(l)):
		print(f'{l[i]:.5g} {vec_units[i]}', end=' ')
	print(units)




