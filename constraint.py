from dataclasses import dataclass


@dataclass
class Constraint:
    position: float  # position along beam [m]
    fix_translation: bool = True  # fix vertical translation
    fix_rotation: bool = False  # fix rotation (slope)

