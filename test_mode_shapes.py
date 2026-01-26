"""Test mode shape computation."""
import numpy as np
from beam_model import BeamModel
from fem import natural_frequencies_and_modes


def test_mode_shapes_returned():
    """Test that mode shapes are computed and returned."""
    model = BeamModel(
        length=2.0,
        radius=0.05,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=7850,
        elements=20,
    )
    
    # Cantilever: fixed at left end (x=0)
    n_nodes = model.elements + 1
    bending_fixed = [0, 1]  # Fix vertical displacement and rotation at node 0
    torsion_fixed = [0]      # Fix torsion at node 0
    
    bend_freq, tors_freq, bend_modes, tors_modes, bf_out, tf_out = natural_frequencies_and_modes(
        model,
        bending_fixed=bending_fixed,
        torsion_fixed=torsion_fixed,
        n_modes=3,
    )
    
    # Verify frequencies are returned
    assert len(bend_freq) == 3
    assert len(tors_freq) == 3
    assert all(f > 0 for f in bend_freq)
    assert all(f > 0 for f in tors_freq)
    
    # Verify mode shapes have correct dimensions
    # Bending: 2 DOFs per node, minus 2 fixed DOFs = 2*21 - 2 = 40 free DOFs
    expected_bend_dofs = 2 * n_nodes - len(bending_fixed)
    assert bend_modes.shape[0] == expected_bend_dofs
    assert bend_modes.shape[1] == 3  # 3 modes
    
    # Torsion: 1 DOF per node, minus 1 fixed DOF = 21 - 1 = 20 free DOFs
    expected_tors_dofs = n_nodes - len(torsion_fixed)
    assert tors_modes.shape[0] == expected_tors_dofs
    assert tors_modes.shape[1] == 3  # 3 modes
    
    # Verify mode shapes are normalized (max value should be 1.0)
    for i in range(3):
        max_bend = np.max(np.abs(bend_modes[:, i]))
        max_tors = np.max(np.abs(tors_modes[:, i]))
        assert np.isclose(max_bend, 1.0, rtol=1e-6)
        assert np.isclose(max_tors, 1.0, rtol=1e-6)
    
    # Verify fixed DOFs are returned
    assert list(bf_out) == bending_fixed
    assert list(tf_out) == torsion_fixed
    
    print("Mode shapes test passed!")
    print(f"Bending frequencies: {bend_freq}")
    print(f"Torsion frequencies: {tors_freq}")
    print(f"Bending mode shape: {bend_modes.shape}")
    print(f"Torsion mode shape: {tors_modes.shape}")


def test_simply_supported_mode_shapes():
    """Test mode shapes for simply supported beam."""
    model = BeamModel(
        length=1.0,
        radius=0.02,
        elastic_modulus=200e9,
        shear_modulus=79.3e9,
        density=7850,
        elements=10,
    )
    
    # Simply supported: fixed vertical displacement at both ends
    n_nodes = model.elements + 1
    bending_fixed = [0, 2 * (n_nodes - 1)]  # Fix vertical at first and last node
    torsion_fixed = [0, n_nodes - 1]        # Fix torsion at both ends
    
    bend_freq, tors_freq, bend_modes, tors_modes, _, _ = natural_frequencies_and_modes(
        model,
        bending_fixed=bending_fixed,
        torsion_fixed=torsion_fixed,
        n_modes=2,
    )
    
    # Verify we get results
    assert len(bend_freq) >= 1
    assert len(tors_freq) >= 1
    assert bend_modes.shape[1] >= 1
    assert tors_modes.shape[1] >= 1
    
    print("Simply supported mode shapes test passed!")


if __name__ == "__main__":
    test_mode_shapes_returned()
    print()
    test_simply_supported_mode_shapes()

