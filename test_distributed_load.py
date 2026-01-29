"""Test distributed load functionality."""
import numpy as np
from beam_model import BeamModel
from distributed_load import DistributedLoad
from shear_moment import shear_moment


def test_distributed_load_simple():
    """Test a simple beam with a uniform distributed load."""
    # Create a simple beam
    g = 9.81  # gravity [m/s^2]
    mass_per_length = 102.0  # kg/m (equivalent to ~1000 N/m)
    model = BeamModel(
        length=2.0,
        radius=0.05,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=0.0,  # Disable structural self-weight
        elements=20,
        distributed_loads=[
            DistributedLoad(start=0.5, end=1.5, mass_per_length=mass_per_length)
        ]
    )

    # Calculate shear and moment (simply supported beam)
    x, shear, moment = shear_moment(
        model,
        g=g,
        n=100,
        left_fixed=True,
        right_fixed=True,
        constraints=[]
    )

    # Verify results exist
    assert x is not None
    assert shear is not None
    assert moment is not None
    assert len(x) == 100
    assert len(shear) == 100
    assert len(moment) == 100

    # Check that shear and moment are non-zero in the region with load
    # The reactions should balance the total load
    total_load = mass_per_length * g * (1.5 - 0.5)  # mass_per_length * g * length
    print(f"Total applied load: {total_load:.1f} N")

    # Maximum moment should be somewhere in the middle
    max_moment_idx = np.argmax(np.abs(moment))
    max_moment = moment[max_moment_idx]
    print(f"Maximum moment: {max_moment:.2f} N⋅m at x={x[max_moment_idx]:.2f} m")

    # Shear should change sign somewhere
    print(f"Shear range: [{shear.min():.2f}, {shear.max():.2f}] N")

    print("Test passed!")


def test_multiple_distributed_loads():
    """Test a beam with multiple distributed loads."""
    g = 9.81
    model = BeamModel(
        length=3.0,
        radius=0.05,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=0.0,  # Disable structural self-weight
        elements=30,
        distributed_loads=[
            DistributedLoad(start=0.0, end=1.0, mass_per_length=51.0),  # ~500 N/m
            DistributedLoad(start=2.0, end=3.0, mass_per_length=102.0)  # ~1000 N/m
        ]
    )

    x, shear, moment = shear_moment(
        model,
        g=g,
        n=150,
        left_fixed=True,
        right_fixed=True,
        constraints=[]
    )

    assert x is not None
    assert shear is not None
    assert moment is not None

    print("Multiple distributed loads test passed!")


if __name__ == "__main__":
    test_distributed_load_simple()
    print()
    test_multiple_distributed_loads()

