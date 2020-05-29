import math
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

class Tractrix:
    def __init__(self, t_0, t_1, p_0, q_0, v, alpha, steps=10000):
        self._t_0 = t_0
        self._t_1 = t_1
        self._steps = steps
        self._dt = (t_1 - t_0) / float(steps)
        self._p_0 = p_0
        self._q_0 = q_0
        self._ell = norm(self._p_0 - self._q_0)
        self._v = v
        self._alpha = alpha

    def _dq_t(self, w_t, t):
        v_t = self._v(t)
        alpha_t = self._alpha(t)
        w_per_t = np.array([- w_t[1], w_t[0]])
        w_norm_t = norm(w_t)

        return (v_t / w_norm_t) * (w_t * math.cos(alpha_t) + w_per_t * math.sin(alpha_t)) * self._dt  

    def evolve(self):  
        w_x = np.zeros(self._steps)
        w_y = np.zeros(self._steps)  
        
        p_x = np.zeros(self._steps)
        p_y = np.zeros(self._steps)
        dp_x = np.zeros(self._steps)
        dp_y = np.zeros(self._steps)
        p_x[0], p_y[0] = self._p_0[0], self._p_0[1]

        q_x = np.zeros(self._steps)
        q_y = np.zeros(self._steps)
        dq_x = np.zeros(self._steps)
        dq_y = np.zeros(self._steps)
        q_x[0], q_y[0] = self._q_0[0], self._q_0[1]

        for i in range(0, self._steps - 1):
            t = i * self._dt
            p_t = np.array([p_x[i], p_y[i]])
            q_t = np.array([q_x[i], q_y[i]])          
            w_t = q_t - p_t
            w_x[i], w_y[i] = w_t[0], w_t[1]            
            dq_t = self._dq_t(w_t, t) 
            dq_x[i], dq_y[i] = dq_t[0], dq_t[1]      
            q_x[i+1], q_y[i+1] = q_t[0] + dq_t[0], q_t[1] + dq_t[1]                        
            dp_t = np.dot(dq_t, w_t) * w_t
            dp_x[i], dp_y[i] = dp_t[0], dp_t[1]   
            p_x[i+1], p_y[i+1] = p_t[0] + dp_t[0], p_t[1] + dp_t[1]

        return p_x, p_y, q_x, q_y, dp_x, dp_y, dq_x, dq_y, w_x, w_y

    def check(self, p_x, p_y, q_x, q_y):
        tolerance = .5
        ell = np.zeros(self._steps - 1)
        scalar = np.zeros(self._steps - 1)
        alpha = np.zeros(self._steps - 1)
        v = np.zeros(self._steps - 1)
        counter = 0
        deviation_too_high = False

        for i in range(0, self._steps - 1):
            t = i * self._dt
            
            _ell = math.sqrt((p_x[i] - q_x[i]) * (p_x[i] - q_x[i]) + (p_y[i] - q_y[i]) * (p_y[i] - q_y[i]))
            ell[i] = _ell
            
            if abs(_ell - self._ell) > tolerance:
                if not deviation_too_high:
                    counter = i
                    deviation_too_high = True
            #assert abs(_ell - self._ell) > tolerance, 'Deviation too high!'

            p_t = np.array([p_x[i], p_y[i]])
            p_s = np.array([p_x[i+1], p_y[i+1]])
            q_t = np.array([q_x[i], q_y[i]]) 
            q_s = np.array([q_x[i+1], q_y[i+1]])
            alpha_t = self._alpha(t)
            v_t = self._v(t)

            dp_t = p_s - p_t 
            dq_t = q_s - q_t
            w_t = q_t - p_t 

            _scalar = dp_t[0] * w_t[1] - dp_t[1] * w_t[0]
            scalar[i] = _scalar

            if abs(_scalar) > tolerance:
                if not deviation_too_high:
                    counter = i
                    deviation_too_high = True
            #assert abs(_scalar) > tolerance, 'Deviation too high!'
            
            _alpha = math.acos(np.dot(w_t, dq_t) / (norm(w_t) * norm(dq_t)))
            alpha[i] = _alpha

            if abs(_alpha - alpha_t) > tolerance:
                if not deviation_too_high:
                    counter = i
                    deviation_too_high = True
            #assert abs(_alpha - alpha_t) > tolerance, 'Deviation too high!'

            _v = norm(dq_t / self._dt)
            v[i] = _v

            if abs(_v - v_t) > tolerance:
                if not deviation_too_high:
                    counter = i
                    deviation_too_high = True
            #assert abs(_v - v_t) > tolerance, 'Deviation too high!'

        return counter, ell, scalar, alpha, v

    def show(self, p_x, p_y, q_x, q_y, dp_x, dp_y, dq_x, dq_y, w_x, w_y):
        for i in range(0, self._steps - 1):
            plt.plot(p_x[0:i], p_y[0:i])
            plt.plot(q_x[0:i], q_y[0:i])

            _w_x = np.linspace(p_x[i], p_x[i] + w_x[i], 100)
            _w_y = np.linspace(p_y[i], p_y[i] + w_y[i], 100)
            plt.plot(_w_x, _w_y)

            _dp_x = np.linspace(p_x[i], p_x[i] + dp_x[i], 100)
            _dp_y = np.linspace(p_y[i], p_y[i] + dp_y[i], 100)
            plt.plot(_dp_x, _dp_y)

            _dq_x = np.linspace(q_x[i], q_x[i] + dq_x[i], 100)
            _dq_y = np.linspace(q_y[i], q_y[i] + dq_y[i], 100)
            plt.plot(_dq_x, _dq_y)

            plt.xlim(-5.0, +5.0)
            plt.ylim(-5.0, +5.0)

            plt.show()

def alpha(t):
    return math.pi / 4.0

def v(t):
    return .25

tractrix = Tractrix(.0, 50.0, np.array([.0, 1.0]), np.array([.0, .0]), v, alpha, 500000)
p_x, p_y, q_x, q_y, dp_x, dp_y, dq_x, dq_y, w_x, w_y = tractrix.evolve()
#tractrix.show(p_x, p_y, q_x, q_y, dp_x, dp_y, dq_x, dq_y, w_x, w_y)

#i, ell, scalar, alpha, v = tractrix.check(p_x, p_y, q_x, q_y)

#plt.plot(ell)
#plt.show()
#plt.plot(scalar)
#plt.show()
#plt.plot(alpha)
#plt.show()
#plt.plot(v)
#plt.show()

plt.plot(p_x, p_y)
plt.plot(q_x, q_y)
plt.axis('equal')
plt.show()
