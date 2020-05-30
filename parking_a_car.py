import math
import numpy as np
import scipy.optimize as op

# p_1 = q_1 + sin(phi)
# p_2 = q_2 + cos(phi)

# p_dot_1 = q_dot_1 + phi_dot * cos(phi)
# p_dot_2 = q_dot_2 - phi_dot * sin(phi)

#     | p_dot_1     sin(phi) | !
# det |                      | = 0
#     | p_dot_2     cos(phi) |

# => 0 = M(q, q_dot, phi, phi_dot) = q_dot_1 * cos(phi) - q_dot_2 * sin(phi) + phi_dot
# 
# L(q, q_dot, phi, phi_dot) = sqrt(2 * q_dot_1 ^ 2 + 2 * q_dot_2 ^ 2 + 2 * phi * (q_dot_1 * cos(phi) - q_dot_2 * sin(phi))) 

def M(q_1, q_2, q_dot_1, q_dot_2, phi, phi_dot):
    return q_dot_1 * math.cos(phi) - q_dot_2 * math.sin(phi) + phi_dot

def L(q_1, q_2, q_dot_1, q_dot_2, phi, phi_dot):
    return math.sqrt(2.0 * q_dot_1 * q_dot_1 + 2.0 * q_dot_2 * q_dot_2 - phi_dot * phi_dot)

class ParkingCar:
    def __init__(self, q_1_0, q_2_0, q_dot_1_0, q_dot_2_0, phi_0, lambda_0, length = 1.0, nbSteps = 5000):
        self._length = length
        self._nbSteps = nbSteps
        self._dt = length / float(nbSteps)
        
        self._q_t = np.zeros((nbSteps, 2))
        self._q_dot_t = np.zeros((nbSteps-1, 2))
        self._q_t[0] = np.array([q_1_0, q_2_0])
        self._q_dot_t[0] = np.array([q_dot_1_0, q_dot_2_0])
        # determine q_t[1] from q_dot[0] and q_t[0]
        self._q_t[1] = self._q_0 + self._dt * self._q_dot_t[0]

        self._phi_t = np.zeros((nbSteps, 2))
        self._phi_dot_t = np.zeros((nbSteps-1, 2))
        self._phi_t[0] = phi_0
        # find phi_dot[0] from 0 = M(q, q_dot, phi, phi_dot)
        self._phi_dot_t[0] = - self._q_dot_1 * math.cos(self._phi_0) + self._q_dot_2 * math.sin(self._phi_0)
        self._phi_t[1] = self._phi_0 + self._dt * self._phi_dot_t[0]

        self._lambda_t = np.zeros((nbSteps-1))
        self._lambda_t[0] = lambda_0

        self._p_t = np.zeros((nbSteps))
        # determine p_t[0] 
        self._p_t[0] = self._q_t[0] + np.array([math.sin(self._phi_0), math.cos(self._phi_0)])

    def _L(self, i):
        return math.sqrt(2.0 * self._q_dot_t[i][0] * self._q_dot_t[i][0] \
            + 2.0 * self._q_dot_t[i][1] * self._q_dot_t[i][1] \
            + self._phi_dot_t[i] * self._phi_dot_t[i] \
            + 2.0 * self._phi_dot_t[i] * (self._q_dot_t[i][0] * math.cos(self._phi_t[i]) \
            - self._q_dot_t[i][1] * math.sin(self._phi_t[i])))

    # dG/dq_dot_1
    def _dG_dq_dot_1(self, i):
        return 2.0 * self._q_dot_t[i][0] / self._L(i) + self._lambda_t[i] * math.cos(self._phi_t[i])

    # dG/dq_dot_2
    def _dG_dq_dot_2(self, i):
        return 2.0 * self._q_dot_t[i][1] / self._L(i) - self._lambda_t[i] * math.sin(self._phi_t[i])

    # dG/dphi
    def _dG_dphi(self, i):
        return - self._lambda_t[i] * (self._q_dot_t[i][0] * math.sin(self._phi_t[i]) \
            - self._q_dot_t[i][1] * math.cos(self._phi_t[i]))

    # dG/dphi_dot
    def _dG_dphi_dot(self, i):
        return - 2.0 * self._phi_dot_t[i] / self._L(i) + self._lambda_t[i]


    def _Lagrange(self, i, q_t_i, lambda_t_i_minus_one):
        # determine q_t_i_minus_one
        self._q_dot_t[i-1] = (q_t_i - self._q_t[i-1]) / self._dt

        # determine phi_t_i_minus_one using M = 0
        self._phi_dot_t[i-1] = - self._q_t[i-1][0] * math.cos(self._phi_t[i-1]) + self._q_t[i-1][1] * math.sin(self._phi_t[i-1])

        # determine phi_t_i
        self._phi_t[i] = self._phi_t[i-1] + self._dt * self._phi_dot_t[i-1]

        # set lambda
        self._lambda_t[i-1] = lambda_t_i_minus_one

        # Lagrange equation for q_1
        lagrange_q_1_i = - (self._dG_dq_dot_1[i-1] - self._dG_dq_dot_1[i-2]) / self._dt

        # Lagrange equation for q_2
        lagrange_q_2_i = - (self._dG_dq_dot_2[i-1] - self._dG_dq_dot_2[i-2]) / self._dt  

        # Lagrange equation for phi
        lagrange_phi_i = self._dG_dphi(i) - (self._dG_dphi_dot(i-1) - self._dG_dphi_dot(i-2)) / self._dt

        return np.array([lagrange_q_1_i, lagrange_q_2_i, lagrange_phi_i])

    def solve(self):
        for t_i in range(2, self._length):
            # initialize start values

            x_0 = [self._q_t[i-1][0], self._q_t[i-1][1], self._phi_t[i-2]]
            op.root(_Lagrange, x0, method='lm', jac=None, tol=None, callback=None, options=None)







