import math 


#grav constant
G = 6.6743 * pow(10,-11)

""" EARTH """

#Sphere of influence in m
SOI_EARTH = 925000000

#Sphere of influence in km
SOI_EARTH_KM = 925000

#mass in kg
MASS_EARTH = 5.974 * pow(10,24)

#radius in m
RADIUS_EARTH = 6378 * pow(10,3)

#radius in km
RADIUS_EARTH_KM = 6378

#gravitational constant mu (m^3/s^2)
MU_EARTH = MASS_EARTH * G

#gravitational constant mu (km^3/s^2)
MU_EARTH_KM = MU_EARTH * pow(10,-9)

EARTH_FLATTENING_FACTOR = (1/298.257)

""" MARS """

#Sphere of influence in m
SOI_MARS = 577000000

#Sphere of influence in km
SOI_MARS_KM = 577000

#mass in kg
MASS_MARS = 641.9 * pow(10,21)

#radius of mars in m 
RADIUS_MARS = 3396 * pow(10,3)

#radius of mars in km
RADIUS_MARS_KM = 3396

#gravitational constant mu (m^3/s^2)
MU_MARS = MASS_MARS * G

#gravitational constant mu (km^3/s^2)
MU_MARS_KM = MU_MARS * pow(10,-9)

""" SUN """

#mass in kg
MASS_SUN = 1.989 * pow(10,30)

#radius in m
RADIUS_SUN = 696000 * pow(10,3)

#radius in km
RADIUS_SUN_KM = 696000

#gravitational constant mu (m^3/s^2)
MU_SUN = MASS_SUN * G

#gravitational constant mu (km^3/s^2)
MU_SUN_KM = MU_SUN * pow(10,-9)

""" JUPITER """

MASS_JUPITER = 1.899 * pow(10,27)

MU_JUPITER = MASS_JUPITER * G

MU_JUPITER_KM = MU_JUPITER * pow(10,-9)

RADIUS_JUPITER_KM = 71490

""" SATURN """

MASS_SATURN = 568.5* pow(10,24)

MU_SATURN = MASS_SATURN * G 

MU_SATURN_KM = MU_SATURN * pow(10, -9)

RADIUS_SATURN_KM = 60270 

""" SEMI-MAJOR AXIS TO SUN """

SUN_EARTH_KM = 149.6 * pow(10,6)

SUN_MARS_KM = 227.9 * pow(10,6)

SUN_JUPITER_KM = 778.6 * pow(10,6)

SUN_SATURN_KM = 1.433 * pow(10, 9)

""" VELOCITIES AROUND SUN (ASSUMES CIRCULAR ORBITS) """

V_EARTH_SUN_KM = math.sqrt(MU_SUN_KM/SUN_EARTH_KM)

V_MARS_SUN_KM = math.sqrt(MU_SUN_KM/SUN_MARS_KM)

V_JUPITER_SUN_KM = math.sqrt(MU_SUN_KM/SUN_JUPITER_KM)

""" SIDE REAL ROTATION PERIOD (s) """

SIDEREAL_ROTATION_MARS = 24.62 * 3600

SIDEREAL_ROTATION_EARTH = 23.9345 * 3600

""" ANGULAR VELOCITIES """

ANG_VEL_EARTH = 72.9217 * pow(10,-6)