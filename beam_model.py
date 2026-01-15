from dataclasses import dataclass, field
from typing import List

import numpy as np

from cross_section import CrossSection
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring


@dataclass
class BeamModel:
    length: float
    radius: float
    elastic_modulus: float
    shear_modulus: float
    density: float
    elements: int = 10
    cross_sections: List[CrossSection] = field(default_factory=list)
    point_masses: List[PointMass] = field(default_factory=list)
    trans_springs: List[TranslationalSpring] = field(default_factory=list)
    tors_springs: List[TorsionalSpring] = field(default_factory=list)

    @property
    def area(self) -> float:
        return self.average_area()

    @property
    def inertia(self) -> float:
        return self.average_inertia()

    @property
    def polar_inertia(self) -> float:
        return self.average_polar_inertia()

    @staticmethod
    def area_for_radius(radius: float) -> float:
        return np.pi * radius**2

    @staticmethod
    def inertia_for_radius(radius: float) -> float:
        return 0.25 * np.pi * radius**4

    @staticmethod
    def polar_inertia_for_radius(radius: float) -> float:
        return 0.5 * np.pi * radius**4

    def radius_at(self, x: float) -> float:
        radius = self.radius
        for section in self.cross_sections:
            if section.start <= x <= section.end:
                radius = section.radius
        return radius

    def element_radii(self) -> np.ndarray:
        if self.elements <= 0 or self.length <= 0:
            return np.array([self.radius])
        L_e = self.length / self.elements
        centers = (np.arange(self.elements) + 0.5) * L_e
        return np.array([self.radius_at(x) for x in centers])

    def average_area(self) -> float:
        if self.length <= 0:
            return self.area_for_radius(self.radius)
        radii = self.element_radii()
        if radii.size == 0:
            return self.area_for_radius(self.radius)
        L_e = self.length / max(self.elements, 1)
        total = sum(self.area_for_radius(r) * L_e for r in radii)
        return total / self.length

    def average_inertia(self) -> float:
        if self.length <= 0:
            return self.inertia_for_radius(self.radius)
        radii = self.element_radii()
        if radii.size == 0:
            return self.inertia_for_radius(self.radius)
        L_e = self.length / max(self.elements, 1)
        total = sum(self.inertia_for_radius(r) * L_e for r in radii)
        return total / self.length

    def average_polar_inertia(self) -> float:
        if self.length <= 0:
            return self.polar_inertia_for_radius(self.radius)
        radii = self.element_radii()
        if radii.size == 0:
            return self.polar_inertia_for_radius(self.radius)
        L_e = self.length / max(self.elements, 1)
        total = sum(self.polar_inertia_for_radius(r) * L_e for r in radii)
        return total / self.length
