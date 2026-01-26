import threading
from typing import Callable, List, Optional, Tuple

import numpy as np

from beam_model import BeamModel
from constraint import Constraint

GRAVITY = 9.81


def _cumulative_trapz(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    if len(x) < 2:
        return np.zeros_like(x)
    dx = np.diff(x)
    avg = 0.5 * (y[1:] + y[:-1])
    return np.concatenate(([0.0], np.cumsum(avg * dx)))


def _distributed_load(model: BeamModel, x: np.ndarray, g: float) -> np.ndarray:
    radii = np.array([model.radius_at(xi) for xi in x])
    areas = np.array([model.area_for_radius(r) for r in radii])
    w = model.density * areas * g

    # Add user-defined distributed loads
    for dist_load in model.distributed_loads:
        # Add load where x is within the range [start, end]
        mask = (x >= dist_load.start) & (x <= dist_load.end)
        w[mask] += dist_load.force_per_length

    return w


ProgressCallback = Callable[[int], None]


def _check_cancel(cancel_event: Optional[threading.Event]) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise RuntimeError("Canceled")


def _report_progress(progress: Optional[ProgressCallback], count: int = 1) -> None:
    if progress is not None:
        progress(count)


def _solve_reactions(
    supports: List[float],
    point_loads: List[Tuple[float, float]],
    x: np.ndarray,
    w: np.ndarray,
    EI: np.ndarray,
    progress: Optional[ProgressCallback] = None,
    cancel_event: Optional[threading.Event] = None,
) -> np.ndarray:
    """Solve for reaction forces at supports using statics and compatibility."""
    n_sup = len(supports)
    if n_sup == 0:
        return np.array([])

    w_total = np.trapezoid(w, x)

    # For one or two supports, use statics (equilibrium equations)
    if n_sup == 1:
        # Single support: reaction equals total load
        total_load = w_total + sum(P for _, P in point_loads)
        _report_progress(progress)
        return np.array([total_load])

    if n_sup == 2:
        # Two supports: use moment and force equilibrium
        x1, x2 = supports[0], supports[1]
        span = x2 - x1

        # Total load and moment about first support
        total_load = w_total
        moment_about_x1 = np.trapezoid(w * (x - x1), x)

        for pos, P in point_loads:
            total_load += P
            moment_about_x1 += P * (pos - x1)

        # R1 + R2 = total_load
        # R2 * span = moment_about_x1
        R2 = moment_about_x1 / span
        R1 = total_load - R2
        _report_progress(progress)
        return np.array([R1, R2])

    # For more than 2 supports: use compatibility (flexibility method)
    # Build flexibility matrix using numerical integration with variable EI.
    F = np.zeros((n_sup, n_sup))
    d = np.zeros(n_sup)

    # Reference: use first and last support as primary supports (simply supported)
    x_left = supports[0]
    x_right = supports[-1]
    span = x_right - x_left
    if span <= 0:
        return np.zeros(n_sup)

    span_mask = (x >= x_left) & (x <= x_right)
    x_span = x[span_mask]
    w_span = w[span_mask]
    EI_span = EI[span_mask]

    if x_span.size < 2:
        return np.zeros(n_sup)

    w_total_span = np.trapezoid(w_span, x_span)
    moment_about_left = np.trapezoid(w_span * (x_span - x_left), x_span)
    for pos, P in point_loads:
        if x_left <= pos <= x_right:
            w_total_span += P
            moment_about_left += P * (pos - x_left)

    R2 = moment_about_left / span
    R1 = w_total_span - R2

    w_cum = _cumulative_trapz(x_span, w_span)
    w_moment_cum = _cumulative_trapz(x_span, w_span * x_span)

    M_load = R1 * (x_span - x_left)
    M_load -= x_span * w_cum - w_moment_cum
    for pos, P in point_loads:
        if x_left <= pos <= x_right:
            M_load -= P * np.maximum(0.0, x_span - pos)

    for i, xi in enumerate(supports):
        _check_cancel(cancel_event)
        if not (x_left <= xi <= x_right):
            continue
        Ra = (x_right - xi) / span
        m_i = np.where(x_span <= xi, Ra * (x_span - x_left), Ra * (x_span - x_left) - (x_span - xi))
        d[i] = np.trapezoid(M_load * m_i / EI_span, x_span)
        _report_progress(progress)
        for j, xj in enumerate(supports):
            _check_cancel(cancel_event)
            if not (x_left <= xj <= x_right):
                continue
            Rb = (x_right - xj) / span
            m_j = np.where(x_span <= xj, Rb * (x_span - x_left), Rb * (x_span - x_left) - (x_span - xj))
            F[i, j] = np.trapezoid(m_i * m_j / EI_span, x_span)
            _report_progress(progress)

    # Solve F * R = -d
    try:
        reactions = np.linalg.solve(F, -d)
    except np.linalg.LinAlgError:
        reactions = np.linalg.lstsq(F, -d, rcond=None)[0]

    _report_progress(progress)
    return reactions


def shear_moment(
    model: BeamModel, g: float = GRAVITY, n: int = 200, left_fixed: bool = True, right_fixed: bool = False,
    constraints: List[Constraint] = None,
    progress: Optional[ProgressCallback] = None,
    cancel_event: Optional[threading.Event] = None,
):
    def _step(count: int = 1) -> None:
        _report_progress(progress, count)

    x = np.linspace(0.0, model.length, n)
    w = _distributed_load(model, x, g)
    radii = np.array([model.radius_at(xi) for xi in x])
    EI = model.elastic_modulus * np.array([model.inertia_for_radius(r) for r in radii])

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
        _step()
        return x, None, None

    supports.sort()
    point_loads = [(pm.position, pm.mass * g) for pm in model.point_masses]
    point_loads.sort(key=lambda p: p[0])

    L = model.length

    # Simple cases first
    if len(supports) == 1:
        # Single support - need to handle specially based on location
        support_pos = supports[0]
        if abs(support_pos) < 1e-9:
            # Left cantilever (fixed at x=0)
            shear = np.zeros_like(x)
            moment = np.zeros_like(x)
            w_cum = _cumulative_trapz(x, w)
            w_moment_cum = _cumulative_trapz(x, w * x)
            for i, xi in enumerate(x):
                _check_cancel(cancel_event)
                load_sum = w_cum[-1] - w_cum[i]
                for pos, load in point_loads:
                    if pos > xi:
                        load_sum += load
                shear[i] = -load_sum

                m = (w_moment_cum[-1] - w_moment_cum[i]) - xi * (w_cum[-1] - w_cum[i])
                for pos, load in point_loads:
                    if pos > xi:
                        m += load * (pos - xi)
                moment[i] = -m
                _step()
            return x, shear, moment
        elif abs(support_pos - L) < 1e-9:
            # Right cantilever (fixed at x=L)
            shear = np.zeros_like(x)
            moment = np.zeros_like(x)
            w_cum = _cumulative_trapz(x, w)
            w_moment_cum = _cumulative_trapz(x, w * x)
            for i, xi in enumerate(x):
                _check_cancel(cancel_event)
                # Loads to the left of xi (same as left cantilever but from left end)
                load_sum = w_cum[i]
                for pos, load in point_loads:
                    if pos < xi:  # Changed from <= to < to match left cantilever logic
                        load_sum += load
                shear[i] = -load_sum  # Negative because loads pull down

                m = xi * w_cum[i] - w_moment_cum[i]
                for pos, load in point_loads:
                    if pos < xi:
                        m += load * (xi - pos)
                moment[i] = -m  # Negative moment
                _step()
            return x, shear, moment
        else:
            # Support at arbitrary position - use general method
            pass  # Fall through to multiple support case

    # Multiple supports: solve for reactions
    reactions = _solve_reactions(
        supports,
        point_loads,
        x,
        w,
        EI,
        progress=progress,
        cancel_event=cancel_event,
    )

    # Calculate shear and moment using equilibrium
    shear = np.zeros_like(x)
    moment = np.zeros_like(x)

    w_cum = _cumulative_trapz(x, w)
    w_moment_cum = _cumulative_trapz(x, w * x)
    for i, xi in enumerate(x):
        _check_cancel(cancel_event)
        # Start with reactions to the left
        shear_val = 0.0
        moment_val = 0.0

        for j, sup_pos in enumerate(supports):
            if sup_pos <= xi:
                shear_val += reactions[j]
                moment_val += reactions[j] * (xi - sup_pos)

        # Subtract distributed load
        shear_val -= w_cum[i]
        moment_val -= xi * w_cum[i] - w_moment_cum[i]

        # Subtract point loads
        for pos, load in point_loads:
            if pos <= xi:
                shear_val -= load
                moment_val -= load * (xi - pos)

        shear[i] = shear_val
        moment[i] = moment_val
        _step()

    return x, shear, moment
