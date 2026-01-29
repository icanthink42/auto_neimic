from dataclasses import dataclass


@dataclass
class DistributedLoad:
    start: float  # Start position along beam [m]
    end: float    # End position along beam [m]
    mass_per_length: float  # Mass per unit length [kg/m]

