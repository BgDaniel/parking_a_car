
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import math



# Problem data.
m = 30
n = 20
np.random.seed(1)
A = np.random.randn(m, n)
b = np.random.randn(m)

# Construct the problem.
x = cp.Variable(n)
objective = cp.Minimize(cp.sum_squares(A*x - b))
constraints = [0 <= x, x <= 1]
prob = cp.Problem(objective, constraints)

# The optimal objective value is returned by `prob.solve()`.
result = prob.solve()
# The optimal value for x is stored in `x.value`.
print(x.value)
# The optimal Lagrange multiplier for a constraint is stored in
# `constraint.dual_value`.
print(constraints[0].dual_value)


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