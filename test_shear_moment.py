import numpy as np

from beam_model import BeamModel
from point_mass import PointMass
from shear_moment import GRAVITY, shear_moment


def _beam():
    return BeamModel(
        length=2.0,
        radius=0.05,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=7850,
        elements=10,
    )


def test_uniform_load_shear_moment_shapes_match_analytical():
    beam = _beam()
    x, V, M = shear_moment(beam, n=5)
    w0 = beam.density * beam.area * GRAVITY
    expected_V = -w0 * (beam.length - x)
    expected_M = -0.5 * w0 * (beam.length - x) ** 2
    np.testing.assert_allclose(V, expected_V, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(M, expected_M, rtol=1e-6, atol=1e-6)


def test_point_mass_load_piecewise():
    beam = _beam()
    beam.density = 0.0  # isolate point load only
    beam.point_masses = [PointMass(position=beam.length / 2, mass=1.0)]
    x, V, M = shear_moment(beam, g=10.0, n=5, left_fixed=True, right_fixed=False)
    P = 10.0
    # expect V = -P for x < mid, 0 for x >= mid
    for xi, v in zip(x, V):
        if xi < beam.length / 2:
            assert np.isclose(v, -P)
        else:
            assert np.isclose(v, 0.0)
    # expect M linear to zero at the load, zero after
    for xi, m in zip(x, M):
        if xi < beam.length / 2:
            expected = -P * (beam.length / 2 - xi)
            assert np.isclose(m, expected)
        else:
            assert np.isclose(m, 0.0)


def test_uniform_load_max_avg():
    beam = _beam()
    x, V, M = shear_moment(beam, n=400, left_fixed=True, right_fixed=False)
    w0 = beam.density * beam.area * GRAVITY
    L = beam.length
    max_s_expected = w0 * L
    avg_s_expected = w0 * L / 2
    max_m_expected = 0.5 * w0 * L**2
    avg_m_expected = w0 * L**2 / 6

    max_s = max(max(V), -min(V))
    avg_s = np.trapezoid(np.abs(V), x) / L
    max_m = max(max(M), -min(M))
    avg_m = np.trapezoid(np.abs(M), x) / L

    np.testing.assert_allclose(max_s, max_s_expected, rtol=1e-3)
    np.testing.assert_allclose(avg_s, avg_s_expected, rtol=1e-3)
    np.testing.assert_allclose(max_m, max_m_expected, rtol=1e-3)
    np.testing.assert_allclose(avg_m, avg_m_expected, rtol=1e-3)


def test_right_fixed_mirrors_left():
    beam = _beam()
    beam.density = 0.0
    beam.point_masses = [PointMass(position=beam.length / 4, mass=1.0)]
    x, Vl, Ml = shear_moment(beam, g=10.0, n=5, left_fixed=True, right_fixed=False)
    x2, Vr, Mr = shear_moment(beam, g=10.0, n=5, left_fixed=False, right_fixed=True)
    assert np.allclose(x, x2)
    # right-fixed: nonzero shear near the right end, zero near the left
    P = 10.0
    assert np.allclose(np.abs(Vr[:2]), 0.0, atol=1e-6)
    assert np.allclose(np.abs(Vr[-2:]), P, atol=1e-6)
    # moments should be largest near the fixed end
    assert np.abs(Mr[-1]) > np.abs(Mr[0])

