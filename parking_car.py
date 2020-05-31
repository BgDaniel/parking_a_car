import numpy as np
import math
from numpy.linalg import inv
import scipy.optimize as op

class ParkingCar:
    def __init__(self, q_start, q_end, phi_start, phi_end, nbSteps = 500):
        self._nbSteps = nbSteps

        self._q_1 = np.zeros((nbSteps))
        self._q_2 = np.zeros((nbSteps))
        self._phi = np.zeros((nbSteps))

        self._q_1[0] = q_start[0]
        self._q_1[-1] = q_end[0]
        self._q_2[0] = q_start[1]
        self._q_2[-1] = q_end[1]
        self._phi[0] = phi_start
        self._phi[-1] = phi_end

        self._q_start = q_start
        self._q_end = q_end
        self._phi_start = phi_start
        self._phi_end = phi_end

    def _phi(self, i):
        return self._phi[i]

    def _dq_1(self, i):
        return self._q_1[i + 1] - self._q_1[i]

    def _dq_2(self, i):
        return self._q_2[i + 1] - self._q_2[i]

    def _dphi(self, i):
        return self._phi[i + 1] - self._phi[i]

    def _dL(self, i):
        return math.sqrt(2.0 * self._dq_1(i) * self._dq_1(i) + 2.0 * self._dq_2(i) * self._dq_2(i)  + self._dphi(i) * self._dphi(i) \
            + 2.0 * self._dphi(i) * (self._dq_1(i) * math.cos(self._phi(i)) - self._dq_2(i) * math.sin(self._phi(i))))

    def _dL_dq_1(self, i):
        if i == 1:
            return - 1.0 / self._dL(1) * (2.0 * self._dq_1(1) + self._dphi(1) * math.cos(self._phi(1))) 
        elif i == self._nbSteps:
            return 1.0 / self._dL(self._nbSteps - 1) * (2.0 * self._dq_1(self._nbSteps - 1) + self._dphi(self._nbSteps - 1) * math.cos(self._phi(self._nbSteps - 1)))
        else:
            return - 1.0 / self._dL(i) * (2.0 * self._dq_1(i) + self._dphi(i) * math.cos(self._phi(i))) \
                + 1.0 / self._dL(i-1) * (2.0 * self._dq_1(i - 1) + self._dphi(i - 1) * math.cos(self._phi(i - 1)))

    def _dL_dq_1(self, i):
        if i == 1:
            return - 1.0 / self._dL(1) * (2.0 * self._dq_2(1) - self._dphi(1) * math.sin(self._phi(1))) 
        elif i == self._nbSteps:
            return 1.0 / self._dL(self._nbSteps - 1) * (2.0 * self._dq_2(self._nbSteps - 1) - self._dphi(self._nbSteps - 1) * math.sin(self._phi(self._nbSteps - 1)))
        else:
            return - 1.0 / self._dL(i) * (2.0 * self._dq_2(i) - self._dphi(i) * math.sin(self._phi(i))) \
                + 1.0 / self._dL(i-1) * (2.0 * self._dq_2(i - 1) - self._dphi(i - 1) * math.sin(self._phi(i - 1)))

    
    def _dL_dphi(self, i):
        if i == 1:
            return 1.0 / (2.0 * self._dL(1)) * (- 2.0 * self._dphi(1) - 2.0 * self._phi(1) * (self._dq_1(1) * math.cos(self._phi(1)) \
                - self._dq_2(1) * math.sin(self._phi(1))) + 2.0 * self._dphi(1) * (- self._dq_1(1) * math.sin(self._phi(1)) - self._dq_2(1) * math.cos(self._phi(1))))
        elif i == N:
            return 1.0 / (2.0 * self._dL(self._nbSteps - 1)) * (2.0 * self._dphi(self._nbSteps - 1) + 2.0 * self._phi(self._nbSteps - 1) * (self._dq_1(self._nbSteps - 1) * math.cos(self._phi(self._nbSteps - 1)) \
                    - self._dq_2(self._nbSteps - 1) * math.sin(self._phi(self._nbSteps - 1))))
        else:
            return 1.0 / (2.0 * self._dL(i)) * (- 2.0 * self._dphi(i) - 2.0 * self._phi(i) * (self._dq_1(i) * math.cos(self._phi(i)) \
                - self._dq_2(i) * math.sin(self._phi(i)) + 2.0 * self._dphi(i) * (- self._dq_1(i) * math.sin(self._phi(i)) - self._dq_2(i) * math.cos(self._phi(i)))) \
                + 1.0 / (2.0 * self._dL(i - 1)) * (2.0 * self._dphi(i - 1) + 2.0 * self._phi(i - 1) * (self._dq_1(i - 1) * math.cos(self._phi(i - 1)) - self._dq_2(i - 1) * math.sin(self._phi(i - 1)))))

    def _dL_dx(self):
        for i in range(1,self._nbSteps):
            dL_dx[i] = self._dL_dq_1(i)
            dL_dx[self._nbSteps + i] = self._dL_dq_2(i)
            dL_dx[2 * self._nbSteps + i] = self._dL_dphi(i)

        return dL_dx

    def _W(self):
        W = np.zeros((3 * self._nbSteps, self._nbSteps - 1))

        for i in range(1, 3 * self._nbSteps):
            for j in range(1, self._nbSteps - 1):
                if 0 < i <= self._nbSteps:
                    if i == j:
                        W[i,j] = - math.cos(self._phi[i])
                    elif i - 1 == j:
                        W[i][j] = math.cos(self._phi[i])
                if self._nbSteps < i <= 2 * self._nbSteps:
                    if i - self._nbSteps == j:
                        W[i][j] = math.sin(self._phi[i - self._nbSteps])
                    elif i - self._nbSteps - 1 == j:
                        W[i][j] = - math.sin(self._phi[i - self._nbSteps])
                if 2 * self._nbSteps < i <= 3 * self._nbSteps:
                    if i - 2 * self._nbSteps == j:
                        W[i][j] = - 1.0 + math.sin(self._phi[i - 2 * self._nbSteps]) * (self._q_1[i - 2 * self._nbSteps+ 1] - self._q_1[i - 2 * self._nbSteps]) \
                            - math.cos(self._phi[i - 2 * self._nbSteps]) * (self._q_2[i - 2 * self._nbSteps + 1] - self._q_2[i - 2 * self._nbSteps])
                    elif i - 2 * self._nbSteps - 1 == j:
                        W[i][j] = 1.0

        return W

    def _projector_W(self):
        W = self._W()
        W_T = np.transpose(W)
        return np.identity(self._n) - W * inv(W_T * W) * W_T

    def _df_projected(self):
        return self._projector_W() * self._dL_dx()

    def park(self):
        q_1_0 = np.linspace(self._q_1[0], self._q_1[-1], num=self._nbSteps)
        q_2_0 = np.linspace(self._q_2[0], self._q_2[-1], num=self._nbSteps)
        phi_0 = np.linspace(self._phi[0], self._phi[-1], num=self._nbSteps)
        x_0 = np.concatenate((q_1_0, q_2_0, phi_0))

        def _F(x):
            for i in range(2, self._nbSteps - 2):
                self._q_1[i] = x[i - 1]
                self._q_2[i] = x[self._nbSteps - 2 + i - 1]
                self._phi[i] = x[2 * self._nbSteps - 4 + i - 1]

            return self._df_projected()

        op.root(_F, x_0, method='lm', jac=None, tol=None, callback=None, options=None)    