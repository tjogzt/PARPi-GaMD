"""
test_gamdSimulation.py

Test the gamdSimulation.py module.
"""

import os
import tempfile
import pytest
import openmm.app as openmm_app
import openmm.unit as unit

from gamd import config
from gamd.gamdSimulation import GamdSimulation, GamdSimulationFactory, load_pdb_positions_and_box_vectors

TEST_DIRECTORY = os.path.dirname(__file__)


class TestGamdSimulation:
    """Test the GamdSimulation class"""

    def test_gamd_simulation_init(self):
        """Test GamdSimulation initialization"""
        sim = GamdSimulation()
        
        # Check that all required attributes are initialized
        assert sim.system is None
        assert sim.integrator is None
        assert sim.simulation is None
        assert sim.traj_reporter is None
        assert sim.first_boost_group is None
        assert sim.second_boost_group is None
        assert sim.first_boost_type is None
        assert sim.second_boost_type is None
        assert sim.topology is None
        assert sim.positions is None
        assert sim.box_vectors is None
        assert sim.platform == "CUDA"
        assert sim.device_index == 0


class TestGamdSimulationFactory:
    """Test the GamdSimulationFactory class"""

    def test_factory_init(self):
        """Test GamdSimulationFactory initialization"""
        factory = GamdSimulationFactory()
        assert factory is not None

    def test_nonbonded_method_validation(self):
        """Test validation of nonbonded methods"""
        factory = GamdSimulationFactory()
        
        # Create a minimal config for testing
        test_config = config.Config()
        
        # Test invalid nonbonded method
        test_config.system.nonbonded_method = "invalid_method"
        
        with pytest.raises(Exception) as excinfo:
            factory.createGamdSimulation(test_config, "CPU", 0)
        
        assert "nonbonded method not found" in str(excinfo.value)

    def test_constraints_validation(self):
        """Test validation of constraints"""
        factory = GamdSimulationFactory()
        
        # Create a minimal config for testing
        test_config = config.Config()
        test_config.system.nonbonded_method = "nocutoff"  # Valid method
        test_config.system.constraints = "invalid_constraint"
        
        with pytest.raises(Exception) as excinfo:
            factory.createGamdSimulation(test_config, "CPU", 0)
        
        assert "constraints not found" in str(excinfo.value)

    def test_algorithm_validation(self):
        """Test validation of integrator algorithms"""
        factory = GamdSimulationFactory()
        
        # Create a minimal config for testing
        test_config = config.Config()
        test_config.system.nonbonded_method = "nocutoff"
        test_config.system.constraints = "none"
        test_config.integrator.algorithm = "invalid_algorithm"
        
        # Mock input files to avoid the "No valid input files found" error
        test_config.input_files.amber = None
        test_config.input_files.charmm = None
        test_config.input_files.gromacs = None
        test_config.input_files.forcefield = None
        
        with pytest.raises(Exception) as excinfo:
            factory.createGamdSimulation(test_config, "CPU", 0)
        
        assert "No valid input files found" in str(excinfo.value)


class TestLoadPdbPositionsAndBoxVectors:
    """Test the load_pdb_positions_and_box_vectors function"""

    def test_load_pdb_function_exists(self):
        """Test that the function exists and is callable"""
        assert callable(load_pdb_positions_and_box_vectors)

    def test_load_pdb_with_invalid_file(self):
        """Test loading PDB with invalid file"""
        with pytest.raises(FileNotFoundError):
            load_pdb_positions_and_box_vectors("nonexistent_file.pdb", False)


@pytest.fixture
def sample_pdb_content():
    """Create a sample PDB file content for testing"""
    pdb_content = """REMARK   1 CREATED WITH OPENMM
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   1      -8.901   4.127  -0.555  1.00  0.00           N
ATOM      2  CA  ALA A   1      -8.608   3.135  -1.618  1.00  0.00           C
ATOM      3  C   ALA A   1      -7.221   2.458  -1.474  1.00  0.00           C
ATOM      4  O   ALA A   1      -6.201   2.698  -2.095  1.00  0.00           O
ATOM      5  CB  ALA A   1      -8.610   3.840  -2.970  1.00  0.00           C
END
"""
    return pdb_content


class TestPDBLoading:
    """Test PDB file loading functionality"""

    def test_load_pdb_with_sample_file(self, sample_pdb_content):
        """Test loading a sample PDB file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(sample_pdb_content)
            f.flush()
            
            try:
                positions, box_vectors = load_pdb_positions_and_box_vectors(f.name, True)
                
                # Check that positions were loaded
                assert positions is not None
                assert hasattr(positions, 'positions')
                
                # Check that box vectors were loaded (from CRYST1 line)
                assert box_vectors is not None
                
            finally:
                os.unlink(f.name)

    def test_load_pdb_without_box_requirement(self, sample_pdb_content):
        """Test loading PDB file without requiring box vectors"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(sample_pdb_content)
            f.flush()
            
            try:
                positions, box_vectors = load_pdb_positions_and_box_vectors(f.name, False)
                
                # Check that positions were loaded
                assert positions is not None
                assert hasattr(positions, 'positions')
                
            finally:
                os.unlink(f.name)


if __name__ == "__main__":
    pytest.main([__file__])
