"""Going to get started with the creation of cascades code."""

import numpy as np
import matplotlib.pyplot as plt


def rk2_custom(f, y0, t, args=()):
    """Solve differential equation using the runge-kutta method.

    Given an equation dy/dt = f(t), f(0) = 0
    https://perso.crans.org/besson/publis/notebooks/Runge-Kutta_methods_for_ODE_integration_in_Python.html

    Args:
        f: Differential equation. Should accept (y, t, *args) IN THAT ORDER.
        y0: Initial value of equation. What does f equal to when t = 0.
        t: Range over which to solve.
        *args: Other arguments for the differential equation.
    Returns:
        Array of values of the solution.
    """
    n = len(t)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(n - 1):
        h = t[i + 1] - t[i]
        k1 = f(y[i], t[i], *args)
        k2 = f(y[i] + (k1 * h) / 2.0, t[i] + (h / 2.0), *args)
        k3 = f(y[i] + (k2 * h) / 2.0, t[i] + (h / 2.0), *args)
        k4 = f(y[i] + (k3 * h), t[i] + h, *args)
        y[i + 1] = y[i] + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


def model1(y, t):
    """dy/dt = 2*t => y(t) = t**2 + C."""
    return 2 * t


initial_value = 100  # Essentially what does y equal to when t is 0?
t_range = np.linspace(-10, 10, 100)
solution = rk2_custom(model1, 100, t_range)

plt.clf()
plt.plot(t_range, solution, label="My Solution")
plt.plot(t_range, [i**2 for i in t_range], "r--", label="Analytical Solution")
plt.title("Runge-Kutta Order 4 Method")
plt.legend()
plt.show()
