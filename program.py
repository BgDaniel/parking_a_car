import numpy as np
from parking_car import ParkingCar

q_start = np.array([1.0, .0])
q_end = np.array([3.0, 2.0])
phi_start = .0
phi_end = .0



n = 101 # points in interval [0,1] to discretize y at
L = 1.5
x = np.arange(n) / float(n)
print(x[1:])
print(x[:-1])

dx = x[1:] - x[:-1]

parkingCar = ParkingCar(q_start, q_end, phi_start, phi_end)
parkingCar.park()

