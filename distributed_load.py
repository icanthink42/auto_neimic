from dataclasses import dataclass


@dataclass
class DistributedLoad:
    start: float  # Start position along beam [m]
    end: float    # End position along beam [m]
    force_per_length: float  # Force per unit length [N/m]

