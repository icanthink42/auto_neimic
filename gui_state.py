from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

from beam_model import BeamModel
from constraint import Constraint
from cross_section import CrossSection
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring


@dataclass
class BeamUIState:
    length: float = 2.0
    radius: float = 0.05
    elastic_modulus: float = 200e9
    shear_modulus: float = 79.3e9
    density: float = 7850
    elements: int = 20
    cross_sections: List[CrossSection] = field(default_factory=list)
    point_masses: List[PointMass] = field(default_factory=list)
    trans_springs: List[TranslationalSpring] = field(default_factory=list)
    tors_springs: List[TorsionalSpring] = field(default_factory=list)
    constraints: List[Constraint] = field(default_factory=list)
    left_fixed: bool = True
    right_fixed: bool = False

    def to_model(self) -> BeamModel:
        return BeamModel(
            length=self.length,
            radius=self.radius,
            elastic_modulus=self.elastic_modulus,
            shear_modulus=self.shear_modulus,
            density=self.density,
            elements=max(4, int(self.elements)),
            cross_sections=list(self.cross_sections),
            point_masses=list(self.point_masses),
            trans_springs=list(self.trans_springs),
            tors_springs=list(self.tors_springs),
        )

    def set_beam(self, params: Tuple[float, float, float, float, float, int]) -> None:
        self.length, self.radius, self.elastic_modulus, self.shear_modulus, self.density, self.elements = params
        self._clamp_cross_sections()

    def _clamp_cross_sections(self) -> None:
        clamped = []
        for section in self.cross_sections:
            start = max(0.0, min(self.length, section.start))
            end = max(0.0, min(self.length, section.end))
            if end <= start:
                continue
            clamped.append(CrossSection(start=start, end=end, radius=section.radius))
        self.cross_sections = clamped
