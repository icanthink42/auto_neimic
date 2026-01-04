import numpy as np

from beam_model import BeamModel
from fem import natural_frequencies
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring


def _cantilever_beam():
    return BeamModel(
        length=2.0,
        radius=0.05,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=7850,
        elements=40,
    )


def test_cantilever_bending_matches_analytical():
    beam = _cantilever_beam()
    bending_bc = [0, 1]
    torsion_bc = [0]
    bend_hz, _ = natural_frequencies(beam, bending_bc, torsion_bc, n_modes=3)

    beta = np.array([1.87510407, 4.69409113, 7.85475744])
    omega = (beta**2) * np.sqrt(
        beam.elastic_modulus * beam.inertia / (beam.density * beam.area * beam.length**4)
    )
    analytical_hz = omega / (2 * np.pi)

    assert np.allclose(bend_hz, analytical_hz, rtol=0.02)


def test_cantilever_torsion_matches_analytical():
    beam = _cantilever_beam()
    bending_bc = [0, 1]
    torsion_bc = [0]
    _, tors_hz = natural_frequencies(beam, bending_bc, torsion_bc, n_modes=3)

    n = np.array([0, 1, 2])
    omega = ((2 * n + 1) * np.pi / (2 * beam.length)) * np.sqrt(
        beam.shear_modulus * beam.polar_inertia / (beam.density * beam.polar_inertia)
    )
    analytical_hz = omega / (2 * np.pi)

    assert np.allclose(tors_hz, analytical_hz, rtol=0.02)


def test_translational_spring_raises_bending_frequency():
    beam = _cantilever_beam()
    bending_bc = [0, 1]
    torsion_bc = [0]
    base_bend, _ = natural_frequencies(beam, bending_bc, torsion_bc, n_modes=1)

    beam.trans_springs = [TranslationalSpring(position=beam.length, k=1e6)]
    with_spring, _ = natural_frequencies(beam, bending_bc, torsion_bc, n_modes=1)

    assert with_spring[0] > base_bend[0]


def test_point_mass_lowers_bending_frequency():
    beam = _cantilever_beam()
    bending_bc = [0, 1]
    torsion_bc = [0]
    base_bend, _ = natural_frequencies(beam, bending_bc, torsion_bc, n_modes=1)

    beam.point_masses = [PointMass(position=beam.length * 0.8, mass=5.0)]
    with_mass, _ = natural_frequencies(beam, bending_bc, torsion_bc, n_modes=1)

    assert with_mass[0] < base_bend[0]

