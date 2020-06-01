import numpy as np
import math
from parking_car import ParkingCar

q_start = np.array([1.0, .0])
q_end = np.array([3.0, 2.0])
phi_start = 5.0 / 4.0 * math.pi
phi_end = 5.0 / 4.0 * math.pi

parkingCar = ParkingCar(q_start, q_end, phi_start, phi_end)
parkingCar.park()

