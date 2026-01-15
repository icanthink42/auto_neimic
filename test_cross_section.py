import numpy as np

from beam_model import BeamModel
from cross_section import CrossSection
from shear_moment import shear_moment


def test_element_radii_respects_cross_sections():
    beam = BeamModel(
        length=1.0,
        radius=0.2,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=7850,
        elements=4,
        cross_sections=[
            CrossSection(start=0.0, end=0.5, radius=0.1),
        ],
    )

    radii = beam.element_radii()
    assert np.allclose(radii, [0.1, 0.1, 0.2, 0.2])


def test_shear_moment_uses_variable_sections():
    beam = BeamModel(
        length=1.0,
        radius=0.2,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=1.0,
        elements=20,
        cross_sections=[
            CrossSection(start=0.0, end=0.5, radius=0.1),
            CrossSection(start=0.5, end=1.0, radius=0.2),
        ],
    )

    x, shear, moment = shear_moment(beam, g=1.0, n=401, left_fixed=True, right_fixed=False)
    area1 = np.pi * 0.1**2
    area2 = np.pi * 0.2**2

    total_load = area1 * 0.5 + area2 * 0.5
    np.testing.assert_allclose(shear[0], -total_load, rtol=2e-3, atol=1e-6)

    target = 0.75
    idx = int(np.argmin(np.abs(x - target)))
    expected_shear = -area2 * (beam.length - x[idx])
    np.testing.assert_allclose(shear[idx], expected_shear, rtol=1e-2, atol=1e-5)

    moment_expected = -(area1 * (0.5**2) / 2 + area2 * (1.0**2 - 0.5**2) / 2)
    np.testing.assert_allclose(moment[0], moment_expected, rtol=1e-2, atol=1e-5)
