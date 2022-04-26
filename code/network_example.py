"""Recreating network example from 2.17 of textbook."""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def model(t, y, rate_constants):
    """Will implement a simple enough solution."""
    a, b, c, d = y
    k1, k2, k3, k4, k5 = rate_constants
    dadt = k1 - k2 * a - k3 * a * b
    dbdt = k2 * a - k3 * a * b
    dcdt = k3 * a * b - k4 * c
    dddt = k3 * a * b - k5 * d
    return [dadt, dbdt, dcdt, dddt]


initial_values = [0, 0, 0, 0]
t = np.linspace(0, 4)
r = [[3, 2, 2.5, 3, 4]]

solution = solve_ivp(model, t_span=[0, 4], t_eval=t, y0=initial_values, args=r)


plt.plot(solution.t, solution.y[0], label="A")
plt.plot(solution.t, solution.y[1], label="B")
plt.plot(solution.t, solution.y[2], label="C")
plt.plot(solution.t, solution.y[3], label="D")
plt.legend()
plt.title("Network Solution")
plt.show()
