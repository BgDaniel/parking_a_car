
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import math

# Problem Data
n = 101 # points in interval [0,1] to discretize y at
L = 1.5

y = cp.Variable(n) # y[i] = y(i/n) for i=0,1,...,n

dy = y[1:] - y[:-1]

x = np.arange(n) / float(n)
dx = x[1:] - x[:-1]



# Construct Optimization Problem
objective = (cp.sum((1.0 + np.multiply(y[1:] - y[:-1], y[1:] - y[:-1]))))


constraints = [
    y[0] == 0, y[-1] == 2
]

prob = cp.Problem(cp.Minimize(objective), constraints)

# Solve Problem
prob.solve(verbose=True)
plt.plot(x, y.value)
plt.show()