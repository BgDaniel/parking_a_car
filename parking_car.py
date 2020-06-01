import numpy as np
import math
from numpy.linalg import inv
import scipy.optimize as op

class ParkingCar:
    def __init__(self, q_start, q_end, phi_start, phi_end, nbSteps = 100):
        self._nbSteps = nbSteps
        self._q_start = q_start
        self._q_end = q_end
        self._phi_start = phi_start
        self._phi_end = phi_end

    def _L(self, q1, q2, phi):
        dq1 = q1[1:] - q1[:-1]
        dq2 = q2[1:] - q2[:-1]
        dphi = phi[1:] - phi[:-1]

        Lq = np.add(np.multiply(dq1, dq1), np.multiply(dq2, dq2))

        H = 2.0 * np.multiply(dphi, np.multiply(dq1, np.cos(phi[:-1])) - np.multiply(dq2, np.sin(phi[:-1])))
        Lp = np.sqrt(np.add(Lq, H))
        Lq = np.sqrt(Lq)

        return Lq.sum() + Lp.sum()

    def _dL(self, q1, q2, phi, eps=10e-6):
        L = self._L(q1, q2, phi)
        
        dLdq1 = np.zeros(self._nbSteps)
        dLdq2 = np.zeros(self._nbSteps)
        dLdphi = np.zeros(self._nbSteps)

        for i in range(0, self._nbSteps):
            delta = np.zeros(self._nbSteps)
            delta[i] = eps
            q1_delta = np.add(q1, delta)
            q2_delta = np.add(q2, delta)
            phi_delta = np.add(phi, delta)

            dLdq1[i] = (self._L(q1_delta, q2, phi) - L) / eps
            dLdq2[i] = (self._L(q1, q2_delta, phi) - L) / eps
            dLdphi[i] = (self._L(q1, q2, phi_delta) - L) / eps

        return np.concatenate((dLdq1, dLdq2, dLdphi))

    def _W(self, q1, q2, phi):
        W = np.zeros(self._nbSteps - 1)
        
        dq1 = q1[1:] - q1[:-1]
        dq2 = q2[1:] - q2[:-1]
        dphi = phi[1:] - phi[:-1]
        cos_phi_q1 = np.multiply(dq1, np.cos(phi[:-1]))
        sin_phi_q2 = np.multiply(dq2, np.sin(phi[:-1]))
        
        return np.add(np.subtract(cos_phi_q1, sin_phi_q2), phi[:-1])

    def _dW(self, q1, q2, phi, eps=10^-5):
        W = self._W(q1, q2, phi)
        dW = np.zeros(self._nbSteps, self._nbSteps - 1)

        for j, W_j in enumerate(W):
            dWjdq1 = np.zeros(self._nbSteps)
            dWjdq2= np.zeros(self._nbSteps)
            dWjdphi = np.zeros(self._nbSteps)
            
            for i in range(0, self._nbSteps):
                delta = np.zeros(self._nbSteps)
                delta[i] = eps
                q1_delta = np.add(q1, delta)
                q2_delta = np.add(q2, delta)
                phi_delta = np.add(phi, delta)

                dWjdq1[i] = (self._L(q1_delta, q2, phi) - W_j) / eps
                dLdq2[i] = (self._L(q1, q2_delta, phi) - W_j) / eps
                dLdphi[i] = (self._L(q1, q2, phi_delta) - W_j) / eps

        return np.concatenate((dLdq1, dLdq2, dLdphi))


        return dW


    def _dW_proj(self, q1, q2, phi):
        dW = self._dW(q1, q2, phi)
        dW_T = np.transpose(dW)
        return np.identity(self._nbSteps) - dW * inv(dW_T * dW) * dW_T

    def _dL_proj(self, q1, q2, phi):
        return self._dW_proj(q1, q2, phi) * self._dL(q1, q2, phi)

    def park(self):
        q1_0 = np.linspace(self._q_start[0], self._q_end[0], num=self._nbSteps)
        q2_0 = np.linspace(self._q_start[1], self._q_end[1], num=self._nbSteps)
        phi_0 = np.linspace(self._phi_start, self._phi_end, num=self._nbSteps)

        dL = self._dL(q1_0, q2_0, phi_0)

        dL_proj = self._dL_proj(q1_0, q2_0, phi_0)

  