import math
import numpy as np

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
    return math.sqrt(2.0 * q_dot_1 * q_dot_1 + 2.0 * q_dot_2 * q_dot_2 + phi_dot * phi_dot /
        + 2.0 * phi_dot * (q_dot_1 * math.cos(phi) - q_dot_2 * math.sin(phi)))

class LagrangeEquation:
    def __init__(self, q_1_0, q_2_0, q_dot_1_0, q_dot_2_0, phi_0, length = 1.0, nbSteps = 5000):
        self._q_0 = np.array([q_1_0, q_2_0])
        self._q_dot_0 = np.array([q_dot_1_0, q_dot_2_0])
        self._phi_0 = phi_0
        self._length = length
        self._nbSteps = nbSteps
        self._dt = length / float(nbSteps)
        
        self._q_t = np.zeros((nbSteps, 2))
        self._q_dot_t = np.zeros((nbSteps-1, 2))
        self._q_t[0] = self._q_0
        self._q_dot_t[0] = self._q_dot_0
        # determine q_t[1] from q_dot[0] and q_t[0]
        self._q_t[1] = self._q_0 + self._dt * self._q_dot_t[0]

        self._phi_t = np.zeros((nbSteps, 2))
        self._phi_dot_t = np.zeros((nbSteps-1, 2))
        self._phi_t[0] = self._phi_0
        # find phi_dot[0] from 0 = M(q, q_dot, phi, phi_dot)
        self._phi_dot_t[0] = - self._q_dot_1 * math.cos(self._phi_0) + self._q_dot_2 * math.sin(self._phi_0)
        self._phi_t[1] = self._phi_0 + self._dt * self._phi_dot_t[0]

        self._p_t = np.zeros((nbSteps))
        # determine p_t[0] 
        self._p_t[0] = self._q_t[0] + np.array([math.sin(self._phi_0), math.cos(self._phi_0)])

    def _L(self, i):
        return L(self._q_t[i][0], self._q_t[i][1], self._q_dot_t[i][0], self._q_dot_t[i][1], self._phi[i], self._phi_dot[i])

    def _M(self, i):
        return M(self._q_t[i][0], self._q_t[i][1], self._q_dot_t[i][0], self._q_dot_t[i][1], self._phi[i], self._phi_dot[i])

    def _dLdq_1_t(self, i):
        return .0

    def _next(self, i, q_dot_i_1, q_dot_i_2, lambda_i_minus_one):
        self._q_dot_t[i] = np.array([q_dot_i_1, q_dot_i_2])
        
        # determine q_t_i_plus_one
        self._q_t[i+1] = self._q_t[i+1] + self._dt * self._q_dot_t[i]

        # determine phi_dot_i
        self._phi_dot_t[i] = - self._q_dot_t[i][0] * math.cos(self._phi_t[i]) + self._q_dot_t[i][1] * math.sin(self._phi_t[i])

        # determine phi_i_plus_one
        self._phi_t[i+1] = self._phi_t[i] + self._dt * self._phi_dot_t_[i]

    def solve(self):
        for t_i in range(1, self._length):
            # solve for q_t_i and lambda_t_i





