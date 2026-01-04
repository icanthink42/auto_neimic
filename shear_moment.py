import numpy as np

from beam_model import BeamModel

GRAVITY = 9.81


def shear_moment(
    model: BeamModel, g: float = GRAVITY, n: int = 200, left_fixed: bool = True, right_fixed: bool = False
):
    x = np.linspace(0.0, model.length, n)
    if (not left_fixed) and (not right_fixed):
        return x, None, None
    w0 = model.density * model.area * g  # distributed load [N/m], downward

    point_loads = [(pm.position, pm.mass * g) for pm in model.point_masses]
    point_loads.sort(key=lambda p: p[0])

    L = model.length

    def _cantilever_left_fixed(xvals, loads):
        shear = np.zeros_like(xvals)
        moment = np.zeros_like(xvals)
        for i, xi in enumerate(xvals):
            load_sum = w0 * (L - xi)
            for pos, load in loads:
                if pos > xi:
                    load_sum += load
            shear[i] = -load_sum

            m = 0.5 * w0 * (L - xi) ** 2
            for pos, load in loads:
                if pos > xi:
                    m += load * (pos - xi)
            moment[i] = -m
        return shear, moment

    if left_fixed and not right_fixed:
        shear, moment = _cantilever_left_fixed(x, point_loads)
        return x, shear, moment

    if right_fixed and not left_fixed:
        # mirror the beam: fixed at right -> compute on reversed coordinates then flip
        x_ref = np.linspace(0.0, L, n)
        loads_ref = [(L - pos, load) for pos, load in point_loads]
        loads_ref.sort(key=lambda p: p[0])
        shear_ref, moment_ref = _cantilever_left_fixed(x_ref, loads_ref)
        shear = -shear_ref[::-1]
        moment = moment_ref[::-1]
        return x, shear, moment

    # both fixed: fixed-end reactions for uniform and point loads
    R1 = w0 * L / 2
    R2 = w0 * L / 2
    M1 = -w0 * L**2 / 12
    M2 = -w0 * L**2 / 12

    for pos, load in point_loads:
        a = pos
        b = L - pos
        R1p = load * (b * b * (3 * a + b)) / L**3
        R2p = load * (a * a * (a + 3 * b)) / L**3
        M1p = -load * a * b * b / (L * L)
        M2p = load * a * a * b / (L * L)
        R1 += R1p
        R2 += R2p
        M1 += M1p
        M2 += M2p

    shear = np.zeros_like(x)
    moment = np.zeros_like(x)
    for i, xi in enumerate(x):
        shear_val = R1 - w0 * xi
        m_val = M1 + R1 * xi - 0.5 * w0 * xi**2
        for pos, load in point_loads:
            if pos <= xi:
                shear_val -= load
                m_val -= load * (xi - pos)
        shear[i] = shear_val
        moment[i] = m_val
    return x, shear, moment

