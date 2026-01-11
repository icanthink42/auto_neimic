import numpy as np
from typing import List, Tuple

from beam_model import BeamModel
from constraint import Constraint

GRAVITY = 9.81


def _solve_reactions(supports: List[float], point_loads: List[Tuple[float, float]],
                     w0: float, L: float, EI: float) -> np.ndarray:
    """Solve for reaction forces at supports using statics and compatibility."""
    n_sup = len(supports)
    if n_sup == 0:
        return np.array([])

    # For one or two supports, use statics (equilibrium equations)
    if n_sup == 1:
        # Single support: reaction equals total load
        total_load = w0 * L + sum(P for _, P in point_loads)
        return np.array([total_load])

    if n_sup == 2:
        # Two supports: use moment and force equilibrium
        x1, x2 = supports[0], supports[1]
        span = x2 - x1

        # Total load and moment about first support
        total_load = w0 * L
        moment_about_x1 = w0 * L * (L / 2 - x1)  # distributed load

        for pos, P in point_loads:
            total_load += P
            moment_about_x1 += P * (pos - x1)

        # R1 + R2 = total_load
        # R2 * span = moment_about_x1
        R2 = moment_about_x1 / span
        R1 = total_load - R2
        return np.array([R1, R2])

    # For more than 2 supports: use compatibility (flexibility method)
    # Build flexibility matrix using simply-supported beam influence coefficients
    F = np.zeros((n_sup, n_sup))
    d = np.zeros(n_sup)

    # Reference: use first and last support as primary supports (simply supported)
    x_left = supports[0]
    x_right = supports[-1]
    span = x_right - x_left

    for i, xi in enumerate(supports):
        # Deflection at xi on a simply supported beam (x_left to x_right) due to loads
        # Position relative to left support
        rel_x = xi - x_left
        a = rel_x / span  # normalized position

        # Deflection due to uniform load on simply supported beam
        if x_left <= xi <= x_right:
            # Uniform load over full span
            w_rel = w0  # load per unit length
            d[i] = -(w_rel * span**4 / (24 * EI)) * a * (1 - a) * (1 + a - 2 * a**2)

        # Deflection due to point loads
        for pos, P in point_loads:
            if x_left <= pos <= x_right:
                b = (pos - x_left) / span
                # Influence coefficient for point load
                if rel_x <= (pos - x_left):
                    delta = (P * span**3 / (6 * EI)) * a * (1 - b) * (2 * b - a - a * b)
                else:
                    delta = (P * span**3 / (6 * EI)) * b * (1 - a) * (2 * a - b - a * b)
                d[i] += delta

        # Flexibility coefficients
        for j, xj in enumerate(supports):
            rel_xj = xj - x_left
            b = rel_xj / span

            if x_left <= xj <= x_right:
                # Influence coefficient for unit load at xj
                if rel_x <= rel_xj:
                    F[i, j] = (span**3 / (6 * EI)) * a * (1 - b) * (2 * b - a - a * b)
                else:
                    F[i, j] = (span**3 / (6 * EI)) * b * (1 - a) * (2 * a - b - a * b)

    # Solve F * R = -d
    try:
        reactions = np.linalg.solve(F, -d)
    except np.linalg.LinAlgError:
        reactions = np.linalg.lstsq(F, -d, rcond=None)[0]

    return reactions


def shear_moment(
    model: BeamModel, g: float = GRAVITY, n: int = 200, left_fixed: bool = True, right_fixed: bool = False,
    constraints: List[Constraint] = None
):
    x = np.linspace(0.0, model.length, n)

    # Collect all support positions that fix translation
    supports = []
    if left_fixed:
        supports.append(0.0)
    if right_fixed:
        supports.append(model.length)
    if constraints:
        for c in constraints:
            if c.fix_translation:
                # Avoid duplicate supports at same location
                if not any(abs(s - c.position) < 1e-9 for s in supports):
                    supports.append(c.position)

    if len(supports) == 0:
        return x, None, None

    supports.sort()
    w0 = model.density * model.area * g  # distributed load [N/m], downward
    point_loads = [(pm.position, pm.mass * g) for pm in model.point_masses]
    point_loads.sort(key=lambda p: p[0])

    L = model.length
    EI = model.elastic_modulus * model.inertia

    # Simple cases first
    if len(supports) == 1:
        # Single support - need to handle specially based on location
        support_pos = supports[0]
        if abs(support_pos) < 1e-9:
            # Left cantilever (fixed at x=0)
            shear = np.zeros_like(x)
            moment = np.zeros_like(x)
            for i, xi in enumerate(x):
                load_sum = w0 * (L - xi)
                for pos, load in point_loads:
                    if pos > xi:
                        load_sum += load
                shear[i] = -load_sum

                m = 0.5 * w0 * (L - xi) ** 2
                for pos, load in point_loads:
                    if pos > xi:
                        m += load * (pos - xi)
                moment[i] = -m
            return x, shear, moment
        elif abs(support_pos - L) < 1e-9:
            # Right cantilever (fixed at x=L)
            shear = np.zeros_like(x)
            moment = np.zeros_like(x)
            for i, xi in enumerate(x):
                # Loads to the left of xi (same as left cantilever but from left end)
                load_sum = w0 * xi
                for pos, load in point_loads:
                    if pos < xi:  # Changed from <= to < to match left cantilever logic
                        load_sum += load
                shear[i] = -load_sum  # Negative because loads pull down

                m = 0.5 * w0 * xi ** 2
                for pos, load in point_loads:
                    if pos < xi:
                        m += load * (xi - pos)
                moment[i] = -m  # Negative moment
            return x, shear, moment
        else:
            # Support at arbitrary position - use general method
            pass  # Fall through to multiple support case

    # Multiple supports: solve for reactions
    reactions = _solve_reactions(supports, point_loads, w0, L, EI)

    # Calculate shear and moment using equilibrium
    shear = np.zeros_like(x)
    moment = np.zeros_like(x)

    for i, xi in enumerate(x):
        # Start with reactions to the left
        shear_val = 0.0
        moment_val = 0.0

        for j, sup_pos in enumerate(supports):
            if sup_pos <= xi:
                shear_val += reactions[j]
                moment_val += reactions[j] * (xi - sup_pos)

        # Subtract distributed load
        shear_val -= w0 * xi
        moment_val -= 0.5 * w0 * xi**2

        # Subtract point loads
        for pos, load in point_loads:
            if pos <= xi:
                shear_val -= load
                moment_val -= load * (xi - pos)

        shear[i] = shear_val
        moment[i] = moment_val

    return x, shear, moment

