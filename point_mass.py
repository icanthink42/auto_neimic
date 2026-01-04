from dataclasses import dataclass


@dataclass
class PointMass:
    position: float
    mass: float
    rotary_inertia: float = 0.0

