
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

# Problem Data
n = 101 # points in interval [0,1] to discretize y at
L = 1.5

y = cp.Variable(n) # y[i] = y(i/n) for i=0,1,...,n

dy = y[1:] - y[:-1]

x = np.arange(n) / float(n)
dx = x[1:] - x[:-1]

dy_dx_sq = np.array(math.sqrt() * dx[i] for i in range(0, len(dx)))
dy_1 = y[1][1:] - y[1][:-1]
dx = x[1:] - x[:-1]

print(dy_0)

# Construct Optimization Problem
objective = (cp.sum((y[:-1] - y[1:]) * ()))


(cp.sum(y[:-1]) + cp.sum(y[1:])) / 2 / n

print(type(objective))

constraints = [
    y[0] == 0, y[-1] == 0,
    cp.sum(cp.norm2(cp.vstack([dy, dx]), axis=0)) <= L,
]

prob = cp.Problem(cp.Maximize(objective), constraints)

# Solve Problem
prob.solve(verbose=True)
plt.plot(x, y.value)
plt.show()