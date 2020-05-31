import nunmpy as np
from numpy.linalg import inv

class DiscreteLagrangian:
    def __init__(self, length = 1.0, nbSteps = 500):
        self._length = length
        self._nbSteps = nbSteps
        self._dt = length / float(nbSteps)
        self._n = 6 * nbSteps - 3
        self._m = 4 * nbSteps - 4 
        self._x = np.zeros((nbSteps))

    def _phi(self, k):
        return self._x[k]

    def _phi_dot(self, k):
        return self._x[self._nbSteps + k]

    def _q_1(self, k):
        return self._x[2 * self._nbSteps - 1 + k]
å‚
    def _q_1_dot(self, k):
        return self._x[3 * self._nbSteps - 1 + k]

    def _q_2(self, k):
        return self._x[4 * self._nbSteps - 2 + k]

    def q_2_dot(self, k):
        return self._x[5 * self._nbSteps - 2 + k]
å‚
    def _df(self, x):
        df = np.zeros((self._n, self._m))

        for i in range(1, self._n):
            if

        return df

    def _W(self, x):
        W = np.zeros((self._n, self._m))

        for i in range(1, self._n):
            for j in range(1, self._m):
                if 0 < i <= self._nbSteps and 0 < j <= self._nbSteps - 1:
                    if i = j:
                        W[i][j] = - 1.0 / self._dt
                    elif i - j == 1:
                        W[i][j] = + 1.0 / self._dt
                elif self._nbSteps < i <= 2 * self._nbSteps - 1 and 0 < j <= self._nbSteps - 1:
                    if i - self._nbSteps == j:
                        W[i][j] = 1.0
                elif 2 * self._nbSteps - 1 < i <= 3 * self._nbSteps - 1 and self._nbSteps < j <= 2 * self._nbSteps - 2:
                    if i - 2 * self._nbSteps + 1 == j - self._nbSteps + 1:
                        W[i][j] = - 1.0 / self._dt
                    elif i - 2 * self._nbSteps == j - self._nbSteps + 1:
                        W[i][j] = + 1.0 / self._dt
                elif 3 * self._nbSteps - 1 < i <= 4 * self._nbSteps - 2:
                    if i - 3 * self._nbSteps + 1 == j - self._nbSteps + 1:
                        W[i][j] = 1.0
                elif 4 * self._nbSteps - 2 < i <= 5 * self._nbSteps - 2 and 2 * self._nbSteps - 2 < j <= 3 * self._nbSteps - 3:
                    if i - 4 * self._nbSteps + 2 == j - 2 * self._nbSteps + 2:
                        W[i][j] = - 1.0 / self._dt
                    elif i - 4 * self._nbSteps + 1 == j - 2 * self._nbSteps + 2:
                        W[i][j] = + 1.0 / self._dt
                elif 5 * self._nbSteps - 2 < i <= self._n and 2 * self._nbSteps - 2 < j <= 3 * self._nbSteps - 3:
                    if i - 5 * self._nbSteps + 2 == j - 2 * self._nbSteps + 2:
                        W[i][j] = 1.0
                elif 0 < i <= self._nbSteps and 3 * self._nbSteps - 3 < j <= self._m:
                    if i == j - 3 * self._nbSteps + 3:
                        W[i][j] = - (self._q_1[i] * math.cos(self._phi[i]) - self._q_2[i] * math.sin(self._phi[i]))
                elif self._nbSteps  < 2 * self._nbSteps - 1 and 3 * self._nbSteps - 3 < j <= self._m:
                    if i - self._nbSteps == j - 3 * self._nbSteps + 3:
                        W[i][j] = 1.0
                elif 3 * self._nbSteps - 1 < i <= 4 * self._nbSteps - 2 and 3 * self._nbSteps - 3 < j <= self._m:
                    if i - 3 * self._nbSteps + 1 = j - 3 * self._nbSteps + 3:
                        W[i][j] = math.cos(self._phi[i - 3 * self._nbSteps + 1])
                elif 5 * self._nbSteps - 2 < i <= self._n and 3 * self._nbSteps - 3 < j <= self._m:
                    if i - 5 * self._nbSteps + 2 = j - 3 * self._nbSteps + 3:
                        W[i][j] = - math.sin(self._phi[i - 5 * self._nbSteps + 2])

        return W

    def projector_W(self, x):
        W = self._W(x)
        W_T = np.transpose(W)
        return np.identity(self._n) - W * inv(W_T * W) * W_T

    def system(self, W):
        return self._projector_W(x) * self._df(x)