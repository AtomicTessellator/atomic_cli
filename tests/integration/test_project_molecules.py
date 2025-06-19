"""
Integration tests for project molecules methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)

To run: uv run pytest tests/integration/test_project_molecules.py -v -m integration
"""

import os
import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.project.molecules import (
    create_project_molecule,
    delete_project_molecule,
    list_project_molecules,
    get_project_molecule,
)


@pytest.fixture(scope="session", autouse=True)
def setup_authentication():
    """Setup authentication for integration tests"""
    load_dotenv()
    username = os.getenv("ATOMICT_USERNAME")
    password = os.getenv("ATOMICT_PASSWORD")
    
    if not username or not password:
        pytest.skip("ATOMICT_USERNAME and ATOMICT_PASSWORD must be set in .env file")
    
    try:
        token = authenticate(username, password)
        os.environ["AT_TOKEN"] = token
        return token
    except Exception as e:
        pytest.skip(f"Authentication failed: {e}")


@pytest.fixture(scope="session")
def test_project_id():
    """Get test project ID from environment"""
    project_id = os.getenv("TEST_PROJECT_ID")
    if not project_id:
        pytest.skip("TEST_PROJECT_ID environment variable must be set")
    return project_id


@pytest.fixture
def cleanup_molecules():
    """Track created molecules for cleanup"""
    created_ids = []
    yield created_ids
    
    # Cleanup after test
    for molecule_id in created_ids:
        try:
            delete_project_molecule(molecule_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def sample_smiles():
    """Sample SMILES strings for testing"""
    return {
        "water": "O",
        "methane": "C",
        "ethanol": "CCO",
        "benzene": "c1ccccc1",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    }


@pytest.mark.integration
class TestProjectMoleculesIntegration:
    """Integration tests for project molecules CRUD operations"""
    
    def test_create_molecule_basic(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test creating a basic project molecule without name"""
        result = create_project_molecule(
            project_id=test_project_id,
            smiles=sample_smiles["water"]
        )
        
        # Track for cleanup
        cleanup_molecules.append(result["id"])
        
        # Verify response structure
        assert "id" in result
        assert result["smiles"] == sample_smiles["water"]
        assert result["project"] == test_project_id
        
    def test_create_molecule_with_name(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test creating a project molecule with name"""
        molecule_name = "Test Water Molecule"
        result = create_project_molecule(
            project_id=test_project_id,
            smiles=sample_smiles["water"],
            name=molecule_name
        )
        
        # Track for cleanup
        cleanup_molecules.append(result["id"])
        
        # Verify response structure
        assert "id" in result
        assert result["smiles"] == sample_smiles["water"]
        assert result["heading"] == molecule_name
        assert result["project"] == test_project_id
        
    def test_create_molecule_complex_smiles(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test creating molecules with complex SMILES strings"""
        result = create_project_molecule(
            project_id=test_project_id,
            smiles=sample_smiles["caffeine"],
            name="Caffeine Test"
        )
        
        # Track for cleanup
        cleanup_molecules.append(result["id"])
        
        # Verify complex SMILES is preserved
        assert result["smiles"] == sample_smiles["caffeine"]
        assert result["heading"] == "Caffeine Test"
        
    def test_get_molecule(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test retrieving a single project molecule"""
        # First create a molecule
        create_result = create_project_molecule(
            project_id=test_project_id,
            smiles=sample_smiles["ethanol"],
            name="Test Ethanol"
        )
        
        # Track for cleanup
        cleanup_molecules.append(create_result["id"])
        molecule_id = create_result["id"]
        
        # Get the molecule
        get_result = get_project_molecule(molecule_id)
        
        # Verify retrieved data matches created data
        assert get_result["id"] == molecule_id
        assert get_result["smiles"] == sample_smiles["ethanol"]
        assert get_result["heading"] == "Test Ethanol"
        assert get_result["project"] == test_project_id
        
    def test_list_molecules_all(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test listing all project molecules"""
        # Create multiple molecules for testing
        molecules_data = [
            {"smiles": sample_smiles["methane"], "name": "Test Methane"},
            {"smiles": sample_smiles["benzene"], "name": "Test Benzene"},
        ]
        
        created_molecules = []
        for mol_data in molecules_data:
            result = create_project_molecule(
                project_id=test_project_id,
                smiles=mol_data["smiles"],
                name=mol_data["name"]
            )
            cleanup_molecules.append(result["id"])
            created_molecules.append(result)
        
        # List all molecules
        list_result = list_project_molecules()
        
        # Verify our created molecules are in the list
        molecule_ids = {mol["id"] for mol in list_result["results"]}
        for created_mol in created_molecules:
            assert created_mol["id"] in molecule_ids
            
    def test_list_molecules_filtered_by_project(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test listing molecules filtered by project"""
        # Create molecules in our test project
        test_molecules = []
        for smiles_name, smiles_val in [("methane", sample_smiles["methane"]), ("water", sample_smiles["water"])]:
            result = create_project_molecule(
                project_id=test_project_id,
                smiles=smiles_val,
                name=f"Test {smiles_name.capitalize()}"
            )
            cleanup_molecules.append(result["id"])
            test_molecules.append(result)
        
        # List molecules filtered by project
        list_result = list_project_molecules(project_id=test_project_id)
        
        # Verify response structure
        assert "results" in list_result
        
        # All returned molecules should belong to our test project
        for mol in list_result["results"]:
            assert mol["project"] == test_project_id
            
        # Our created molecules should be in the filtered list
        filtered_molecule_ids = {mol["id"] for mol in list_result["results"]}
        for test_mol in test_molecules:
            assert test_mol["id"] in filtered_molecule_ids
    
    def test_delete_molecule(self, test_project_id, sample_smiles):
        """Test deleting a project molecule"""
        # First create a molecule
        create_result = create_project_molecule(
            project_id=test_project_id,
            smiles=sample_smiles["methane"],
            name="To Delete"
        )
        
        molecule_id = create_result["id"]
        
        # Verify it exists
        get_result = get_project_molecule(molecule_id)
        assert get_result["id"] == molecule_id
        
        # Delete it
        delete_result = delete_project_molecule(molecule_id)
        # Delete endpoint may return empty response or success message
        
        # Verify it's deleted (should raise exception or return 404)
        with pytest.raises(Exception):
            get_project_molecule(molecule_id)


@pytest.mark.integration
class TestProjectMoleculesWorkflow:
    """End-to-end workflow tests for project molecules"""
    
    def test_complete_molecule_crud_workflow(self, test_project_id, sample_smiles):
        """Test complete molecule CRUD workflow"""
        molecule_name = "Workflow Test Molecule"
        
        # Step 1: Create molecule
        create_result = create_project_molecule(
            project_id=test_project_id,
            smiles=sample_smiles["benzene"],
            name=molecule_name
        )
        
        molecule_id = create_result["id"]
        assert create_result["smiles"] == sample_smiles["benzene"]
        assert create_result["heading"] == molecule_name
        
        # Step 2: Get molecule details
        get_result = get_project_molecule(molecule_id)
        assert get_result["id"] == molecule_id
        assert get_result["heading"] == molecule_name
        assert get_result["smiles"] == sample_smiles["benzene"]
        
        # Step 3: Verify it appears in project molecule list
        list_result = list_project_molecules(project_id=test_project_id)
        molecule_ids = {mol["id"] for mol in list_result["results"]}
        assert molecule_id in molecule_ids
        
        # Step 4: Verify it appears in all molecules list
        all_list_result = list_project_molecules()
        all_molecule_ids = {mol["id"] for mol in all_list_result["results"]}
        assert molecule_id in all_molecule_ids
        
        # Step 5: Delete molecule
        delete_result = delete_project_molecule(molecule_id)
        
        # Step 6: Verify deletion (should raise exception or return 404)
        with pytest.raises(Exception):
            get_project_molecule(molecule_id)
            
    def test_multiple_molecules_same_project(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test creating multiple molecules in the same project"""
        molecules_to_create = [
            {"smiles": sample_smiles["water"], "name": "Water Test"},
            {"smiles": sample_smiles["methane"], "name": "Methane Test"},
            {"smiles": sample_smiles["ethanol"], "name": "Ethanol Test"},
        ]
        
        created_molecules = []
        
        # Create all molecules
        for mol_data in molecules_to_create:
            result = create_project_molecule(
                project_id=test_project_id,
                smiles=mol_data["smiles"],
                name=mol_data["name"]
            )
            cleanup_molecules.append(result["id"])
            created_molecules.append(result)
        
        # Verify all molecules exist individually
        for created_mol in created_molecules:
            get_result = get_project_molecule(created_mol["id"])
            assert get_result["id"] == created_mol["id"]
            assert get_result["project"] == test_project_id
        
        # Verify all molecules appear in project list
        list_result = list_project_molecules(project_id=test_project_id)
        listed_molecule_ids = {mol["id"] for mol in list_result["results"]}
        
        for created_mol in created_molecules:
            assert created_mol["id"] in listed_molecule_ids
            
    def test_molecule_smiles_preservation(self, test_project_id, sample_smiles, cleanup_molecules):
        """Test that various SMILES formats are preserved correctly"""
        test_cases = [
            {"name": "Simple", "smiles": sample_smiles["water"]},
            {"name": "Aromatic", "smiles": sample_smiles["benzene"]},
            {"name": "Complex", "smiles": sample_smiles["caffeine"]},
        ]
        
        for test_case in test_cases:
            # Create molecule
            create_result = create_project_molecule(
                project_id=test_project_id,
                smiles=test_case["smiles"],
                name=test_case["name"]
            )
            
            cleanup_molecules.append(create_result["id"])
            
            # Verify SMILES is preserved exactly
            assert create_result["smiles"] == test_case["smiles"]
            
            # Verify through get endpoint as well
            get_result = get_project_molecule(create_result["id"])
            assert get_result["smiles"] == test_case["smiles"]
