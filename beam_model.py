from dataclasses import dataclass, field
from typing import List

import numpy as np

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
    point_masses: List[PointMass] = field(default_factory=list)
    trans_springs: List[TranslationalSpring] = field(default_factory=list)
    tors_springs: List[TorsionalSpring] = field(default_factory=list)

    @property
    def area(self) -> float:
        return np.pi * self.radius**2

    @property
    def inertia(self) -> float:
        return 0.25 * np.pi * self.radius**4

    @property
    def polar_inertia(self) -> float:
        return 0.5 * np.pi * self.radius**4

