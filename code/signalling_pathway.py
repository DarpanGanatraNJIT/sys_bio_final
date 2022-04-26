"""Creating the signalling pathway from Figure 6.2."""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def signal_pathway_model(t, y, rate_constants):
    """Create signalling pathway from Figure 6.2.

    =============================
    R + L <--[k1/k11]-> RL
    P + RL ----[k2]---> P* + RL
    P* --------[k3]---> P
    =============================
    R: Receptor
    L: Ligand
    RL: Active receptor-ligand complex
    P, P*: Inactive and active response regulator proteins
    """
    r, rl, p, ps = y
    k1, k2, k3, k11 = rate_constants
    # Introduce the concentration of ligand here
    if t >= 1 and t <= 3:
        l = 3
    else:
        l = 0
    drdt = -k1 * r * l + k11 * rl
    drldt = k1 * r * l - k11 * rl
    dpdt = -k2 * p * rl + k3 * ps
    dpsdt = k2 * p * rl - k3 * ps
    return [drdt, drldt, dpdt, dpsdt]


initial_values = [2, 0, 8, 0]
rate_constants = [[5, 6, 2, 1]]
t_values = np.linspace(0, 10, num=100)

ligand_concentration = np.zeros(len(t_values))
ligand_concentration[np.where((t_values >= 1) & (t_values <= 3))] = 3


solution = solve_ivp(
    signal_pathway_model,
    t_span=[0, 10],
    t_eval=t_values,
    y0=initial_values,
    args=rate_constants,
    dense_output=True,
    rtol=1e-10,
    atol=1e-10,
)

plt.clf()
plt.plot(solution.t, solution.y[3], label="Active Protein (P*)")
plt.plot(solution.t, solution.y[1], "-.", label="Receptor-Ligand Complex (RL)")
plt.plot(solution.t, ligand_concentration, ".", label="Total Ligand (L)")
plt.legend()
plt.show()
