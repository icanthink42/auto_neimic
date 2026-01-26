from __future__ import annotations

from typing import Callable, Optional, Sequence, Tuple

import threading

import numpy as np

from beam_model import BeamModel
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring


def beam_nodes(model: BeamModel) -> np.ndarray:
    return np.linspace(0.0, model.length, model.elements + 1)


def _shape_factors(x: float, x1: float, x2: float) -> Tuple[float, float]:
    xi = (x - x1) / (x2 - x1)
    xi = np.clip(xi, 0.0, 1.0)
    return 1.0 - xi, xi


ProgressCallback = Callable[[int], None]


def _check_cancel(cancel_event: Optional[threading.Event]) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise RuntimeError("Canceled")


def _report_progress(progress: Optional[ProgressCallback], count: int = 1) -> None:
    if progress is not None:
        progress(count)


def assemble_bending(
    model: BeamModel,
    progress: Optional[ProgressCallback] = None,
    cancel_event: Optional[threading.Event] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    n_nodes = model.elements + 1
    dof = 2 * n_nodes
    K = np.zeros((dof, dof))
    M = np.zeros((dof, dof))

    L_e = model.length / model.elements
    k_template = np.array(
        [
            [12, 6 * L_e, -12, 6 * L_e],
            [6 * L_e, 4 * L_e**2, -6 * L_e, 2 * L_e**2],
            [-12, -6 * L_e, 12, -6 * L_e],
            [6 * L_e, 2 * L_e**2, -6 * L_e, 4 * L_e**2],
        ]
    )
    m_template = np.array(
        [
            [156, 22 * L_e, 54, -13 * L_e],
            [22 * L_e, 4 * L_e**2, 13 * L_e, -3 * L_e**2],
            [54, 13 * L_e, 156, -22 * L_e],
            [-13 * L_e, -3 * L_e**2, -22 * L_e, 4 * L_e**2],
        ]
    )

    radii = model.element_radii()
    for e in range(model.elements):
        _check_cancel(cancel_event)
        radius = radii[e] if e < len(radii) else model.radius
        EI = model.elastic_modulus * model.inertia_for_radius(radius)
        rhoA = model.density * model.area_for_radius(radius)
        k_local = (EI / L_e**3) * k_template
        m_local = (rhoA * L_e / 420) * m_template
        idx = 2 * e
        dofs = [idx, idx + 1, idx + 2, idx + 3]
        K[np.ix_(dofs, dofs)] += k_local
        M[np.ix_(dofs, dofs)] += m_local
        _report_progress(progress)

    x_nodes = beam_nodes(model)
    for pm in model.point_masses:
        _check_cancel(cancel_event)
        e = min(np.searchsorted(x_nodes, pm.position) - 1, model.elements - 1)
        e = max(e, 0)
        x1, x2 = x_nodes[e], x_nodes[e + 1]
        N1, N2 = _shape_factors(pm.position, x1, x2)
        add = pm.mass * np.array([[N1**2, N1 * N2], [N1 * N2, N2**2]])
        w_dofs = [2 * e, 2 * (e + 1)]
        M[np.ix_(w_dofs, w_dofs)] += add
        _report_progress(progress)

    for sp in model.trans_springs:
        _check_cancel(cancel_event)
        e = min(np.searchsorted(x_nodes, sp.position) - 1, model.elements - 1)
        e = max(e, 0)
        x1, x2 = x_nodes[e], x_nodes[e + 1]
        N1, N2 = _shape_factors(sp.position, x1, x2)
        add = sp.k * np.array([[N1**2, N1 * N2], [N1 * N2, N2**2]])
        w_dofs = [2 * e, 2 * (e + 1)]
        K[np.ix_(w_dofs, w_dofs)] += add
        _report_progress(progress)

    return K, M


def assemble_torsion(
    model: BeamModel,
    progress: Optional[ProgressCallback] = None,
    cancel_event: Optional[threading.Event] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    n_nodes = model.elements + 1
    K = np.zeros((n_nodes, n_nodes))
    M = np.zeros((n_nodes, n_nodes))

    L_e = model.length / model.elements
    k_template = np.array([[1, -1], [-1, 1]])
    m_template = np.array([[2, 1], [1, 2]])

    radii = model.element_radii()
    for e in range(model.elements):
        _check_cancel(cancel_event)
        radius = radii[e] if e < len(radii) else model.radius
        GJ = model.shear_modulus * model.polar_inertia_for_radius(radius)
        rhoJ = model.density * model.polar_inertia_for_radius(radius)
        k_local = (GJ / L_e) * k_template
        m_local = (rhoJ * L_e / 6) * m_template
        dofs = [e, e + 1]
        K[np.ix_(dofs, dofs)] += k_local
        M[np.ix_(dofs, dofs)] += m_local
        _report_progress(progress)

    x_nodes = beam_nodes(model)
    for pm in model.point_masses:
        _check_cancel(cancel_event)
        if pm.rotary_inertia == 0:
            _report_progress(progress)
            continue
        e = min(np.searchsorted(x_nodes, pm.position) - 1, model.elements - 1)
        e = max(e, 0)
        x1, x2 = x_nodes[e], x_nodes[e + 1]
        N1, N2 = _shape_factors(pm.position, x1, x2)
        add = pm.rotary_inertia * np.array([[N1**2, N1 * N2], [N1 * N2, N2**2]])
        dofs = [e, e + 1]
        M[np.ix_(dofs, dofs)] += add
        _report_progress(progress)

    for sp in model.tors_springs:
        _check_cancel(cancel_event)
        e = min(np.searchsorted(x_nodes, sp.position) - 1, model.elements - 1)
        e = max(e, 0)
        x1, x2 = x_nodes[e], x_nodes[e + 1]
        N1, N2 = _shape_factors(sp.position, x1, x2)
        add = sp.k * np.array([[N1**2, N1 * N2], [N1 * N2, N2**2]])
        dofs = [e, e + 1]
        K[np.ix_(dofs, dofs)] += add
        _report_progress(progress)

    return K, M


def apply_boundary_conditions(
    K: np.ndarray, M: np.ndarray, fixed_dofs: Sequence[int]
) -> Tuple[np.ndarray, np.ndarray]:
    keep = np.array([i for i in range(K.shape[0]) if i not in fixed_dofs])
    return K[np.ix_(keep, keep)], M[np.ix_(keep, keep)]


def solve_frequencies(K: np.ndarray, M: np.ndarray, n_modes: int) -> np.ndarray:
    A = np.linalg.solve(M, K)
    eigvals, _ = np.linalg.eig(A)
    eigvals = np.real(eigvals[eigvals > 0])
    eigvals.sort()
    omegas = np.sqrt(eigvals)
    return omegas[: min(n_modes, len(omegas))] / (2 * np.pi)


def solve_modes(
    K: np.ndarray, M: np.ndarray, n_modes: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve for natural frequencies and mode shapes.
    
    Returns:
        frequencies: Natural frequencies in Hz
        mode_shapes: Mode shape vectors (columns are modes)
    """
    A = np.linalg.solve(M, K)
    eigvals, eigvecs = np.linalg.eig(A)
    
    # Filter positive eigenvalues
    mask = eigvals > 0
    eigvals = np.real(eigvals[mask])
    eigvecs = np.real(eigvecs[:, mask])
    
    # Sort by eigenvalue
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    # Convert to frequencies
    omegas = np.sqrt(eigvals) / (2 * np.pi)
    
    # Normalize mode shapes by maximum displacement
    for i in range(eigvecs.shape[1]):
        max_val = np.max(np.abs(eigvecs[:, i]))
        if max_val > 1e-10:
            eigvecs[:, i] /= max_val
    
    n = min(n_modes, len(omegas))
    return omegas[:n], eigvecs[:, :n]


def natural_frequencies(
    model: BeamModel,
    bending_fixed: Sequence[int],
    torsion_fixed: Sequence[int],
    n_modes: int = 6,
    progress: Optional[ProgressCallback] = None,
    cancel_event: Optional[threading.Event] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    Kb, Mb = assemble_bending(model, progress=progress, cancel_event=cancel_event)
    Kt, Mt = assemble_torsion(model, progress=progress, cancel_event=cancel_event)

    Kb_r, Mb_r = apply_boundary_conditions(Kb, Mb, bending_fixed)
    Kt_r, Mt_r = apply_boundary_conditions(Kt, Mt, torsion_fixed)

    bend = solve_frequencies(Kb_r, Mb_r, n_modes)
    _report_progress(progress)
    tors = solve_frequencies(Kt_r, Mt_r, n_modes)
    _report_progress(progress)
    return bend, tors


def natural_frequencies_and_modes(
    model: BeamModel,
    bending_fixed: Sequence[int],
    torsion_fixed: Sequence[int],
    n_modes: int = 6,
    progress: Optional[ProgressCallback] = None,
    cancel_event: Optional[threading.Event] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Sequence[int], Sequence[int]]:
    """Compute natural frequencies and mode shapes.
    
    Returns:
        bend_freq: Bending natural frequencies [Hz]
        tors_freq: Torsion natural frequencies [Hz]
        bend_modes: Bending mode shapes (reduced DOFs)
        tors_modes: Torsion mode shapes (reduced DOFs)
        bending_fixed: Fixed bending DOFs (for reconstruction)
        torsion_fixed: Fixed torsion DOFs (for reconstruction)
    """
    Kb, Mb = assemble_bending(model, progress=progress, cancel_event=cancel_event)
    Kt, Mt = assemble_torsion(model, progress=progress, cancel_event=cancel_event)

    Kb_r, Mb_r = apply_boundary_conditions(Kb, Mb, bending_fixed)
    Kt_r, Mt_r = apply_boundary_conditions(Kt, Mt, torsion_fixed)

    bend_freq, bend_modes = solve_modes(Kb_r, Mb_r, n_modes)
    _report_progress(progress)
    tors_freq, tors_modes = solve_modes(Kt_r, Mt_r, n_modes)
    _report_progress(progress)
    
    return bend_freq, tors_freq, bend_modes, tors_modes, bending_fixed, torsion_fixed
