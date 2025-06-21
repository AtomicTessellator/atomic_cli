"""
Integration tests for HT ML Relaxation that includes creating the prerequisite HT SQS exploration.

This test creates its own HT SQS exploration first, then tests HT MLRelax functionality.

These tests require:
- Valid authentication credentials (ATOMICT_USERNAME, ATOMICT_PASSWORD)
- A test project (TEST_PROJECT_ID)
- A test structure (TEST_STRUCTURE_ID, TEST_STRUCTURE_TYPE)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_ht_mlrelax_with_setup.py -v -m integration
"""

import os
from typing import List

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate

from atomict.simulation.ht_sqs import (
    create_ht_sqs_exploration,
    delete_ht_sqs_exploration,
    get_ht_sqs_exploration,
)

from atomict.simulation.ht_mlrelax import (
    ML_MODELS,
    create_ht_mlrelax_exploration,
    delete_ht_mlrelax_exploration,
    get_ht_mlrelax_exploration,
    list_ht_mlrelax_explorations,
    update_ht_mlrelax_exploration,
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
    """Test project ID from environment"""
    project_id = os.getenv("TEST_PROJECT_ID")
    if not project_id:
        pytest.skip("TEST_PROJECT_ID environment variable not set")
    return project_id


@pytest.fixture(scope="session")
def test_structure_id():
    """Test structure ID from environment"""
    structure_id = os.getenv("TEST_STRUCTURE_ID")
    if not structure_id:
        pytest.skip("TEST_STRUCTURE_ID environment variable not set")
    return structure_id


@pytest.fixture(scope="session")
def test_structure_type():
    """Test structure type from environment (defaults to userupload)"""
    return os.getenv("TEST_STRUCTURE_TYPE", "userupload")


@pytest.fixture(scope="session")
def ht_sqs_exploration_for_mlrelax(test_project_id, test_structure_id, test_structure_type):
    """Create an HT SQS exploration specifically for testing HT MLRelax"""
    
    # Create HT SQS exploration
    ht_sqs_result = create_ht_sqs_exploration(
        project_id=test_project_id,
        name="HT SQS for MLRelax Integration Test",
        target_concentrations=[
            {"element": "Al", "weight": 0.7},
            {"element": "Ni", "weight": 0.3}
        ],
        generated_permutations=[
            {"Al": 0.8, "Ni": 0.2},
            {"Al": 0.6, "Ni": 0.4}
        ],
        structure_id=test_structure_id,
        structure_type=test_structure_type,
        num_structures_per_permutation=1,
        description="Created for HT MLRelax integration testing"
    )
    
    print(f"Created HT SQS exploration: {ht_sqs_result}")
    
    # Find the created exploration by listing and matching the name
    from atomict.simulation.ht_sqs import list_ht_sqs_explorations
    
    list_result = list_ht_sqs_explorations(project_id=test_project_id, limit=50)
    ht_sqs_exploration = None
    
    for exp in list_result["results"]:
        if exp["name"] == "HT SQS for MLRelax Integration Test":
            ht_sqs_exploration = exp
            break
    
    if not ht_sqs_exploration:
        pytest.fail("Could not find created HT SQS exploration")
    
    ht_sqs_id = ht_sqs_exploration["id"]
    print(f"Found HT SQS exploration ID: {ht_sqs_id}")
    
    yield ht_sqs_id
    
    # Cleanup - delete the HT SQS exploration
    try:
        delete_ht_sqs_exploration(ht_sqs_id)
        print(f"Cleaned up HT SQS exploration: {ht_sqs_id}")
    except Exception as e:
        print(f"Failed to cleanup HT SQS exploration {ht_sqs_id}: {e}")


@pytest.fixture
def cleanup_ht_mlrelax_explorations():
    """Track created HT ML Relaxation explorations for cleanup"""
    created_ids = []

    def track_id(exploration_id):
        created_ids.append(exploration_id)

    yield track_id

    # Cleanup
    for exploration_id in created_ids:
        try:
            delete_ht_mlrelax_exploration(exploration_id)
            print(f"Cleaned up HT MLRelax exploration: {exploration_id}")
        except Exception as e:
            print(f"Failed to cleanup HT MLRelax exploration {exploration_id}: {e}")


@pytest.mark.integration
class TestHTMLRelaxIntegrationWithSetup:
    """Integration tests for HT ML Relaxation with complete setup"""

    def test_create_ht_mlrelax_exploration_basic(
        self,
        test_project_id,
        ht_sqs_exploration_for_mlrelax,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test creating a HT ML Relaxation exploration with basic parameters"""

        result = create_ht_mlrelax_exploration(
            project_id=test_project_id,
            name="Integration Test HT MLRelax Basic",
            source_ht_sqs_exploration_id=ht_sqs_exploration_for_mlrelax,
            f_max=0.05,
            model="mattersim_1_0_0_5m",
        )

        print(f"Create result: {result}")

        # Should return success message with count
        assert "message" in result
        assert "relaxation" in result["message"].lower()

    def test_create_ht_mlrelax_exploration_all_models(
        self,
        test_project_id,
        ht_sqs_exploration_for_mlrelax,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test creating HT ML Relaxation explorations with all ML models"""

        for model_name in ML_MODELS.keys():
            result = create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name=f"Integration Test HT MLRelax {model_name}",
                source_ht_sqs_exploration_id=ht_sqs_exploration_for_mlrelax,
                model=model_name,
                f_max=0.03,
                description=f"Testing {model_name} model",
            )

            print(f"Model {model_name} result: {result}")
            
            assert "message" in result
            assert "relaxation" in result["message"].lower()

    def test_ht_mlrelax_full_workflow(
        self,
        test_project_id,
        ht_sqs_exploration_for_mlrelax,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test complete HT ML Relaxation workflow: create -> list -> get -> update"""

        # 1. Create HT ML Relaxation exploration
        create_result = create_ht_mlrelax_exploration(
            project_id=test_project_id,
            name="Integration Test HT MLRelax Full Workflow",
            description="Full workflow integration test",
            source_ht_sqs_exploration_id=ht_sqs_exploration_for_mlrelax,
            f_max=0.04,
            model="orb_v3_conservative",
        )

        print(f"Create workflow result: {create_result}")
        assert "relaxation" in create_result["message"].lower()

        # 2. List explorations and find our creation
        list_result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=50)
        print(f"Found {len(list_result['results'])} HT MLRelax explorations")

        # Find our exploration by name
        our_exploration = None
        for exp in list_result["results"]:
            if exp["name"] == "Integration Test HT MLRelax Full Workflow":
                our_exploration = exp
                break

        assert our_exploration is not None, "Created exploration not found in list"
        cleanup_ht_mlrelax_explorations(our_exploration["id"])

        print(f"Found our exploration: {our_exploration['id']}")

        # 3. Get exploration details
        detail_result = get_ht_mlrelax_exploration(our_exploration["id"])

        assert "id" in detail_result
        assert detail_result["name"] == "Integration Test HT MLRelax Full Workflow"
        assert detail_result["description"] == "Full workflow integration test"
        assert detail_result["f_max"] == 0.04
        assert detail_result["model"] == ML_MODELS["orb_v3_conservative"]

        print(f"Detailed exploration: {detail_result}")

        # 4. Get with children to see individual ML relaxations
        detail_with_children = get_ht_mlrelax_exploration(
            our_exploration["id"], include_children=True
        )

        print(f"Children count: {len(detail_with_children.get('children', []))}")
        
        if detail_with_children.get("children"):
            child = detail_with_children["children"][0]
            print(f"First child: {child}")
            assert "mlrelax" in child
            assert child["mlrelax"] is not None

        # 5. Update exploration parameters
        update_result = update_ht_mlrelax_exploration(
            our_exploration["id"],
            name="Updated Integration Test Name",
            f_max=0.02,
            model="mattersim_1_0_0_5m",
            description="Updated description",
        )

        print(f"Update result: {update_result}")

        assert update_result["name"] == "Updated Integration Test Name"
        assert update_result["f_max"] == 0.02
        assert update_result["model"] == ML_MODELS["mattersim_1_0_0_5m"]
        assert update_result["description"] == "Updated description"

    def test_ht_mlrelax_validation(
        self,
        test_project_id,
        ht_sqs_exploration_for_mlrelax,
    ):
        """Test validation and error handling"""

        # Test invalid model validation (should fail before API call)
        with pytest.raises(ValueError, match="Invalid model"):
            create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name="Invalid Model Test",
                source_ht_sqs_exploration_id=ht_sqs_exploration_for_mlrelax,
                model="invalid_model_name",
            )

        # Test invalid exploration ID (should fail at API level)
        with pytest.raises(Exception):  # Will be API error
            create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name="Invalid SQS ID Test",
                source_ht_sqs_exploration_id="00000000-0000-0000-0000-000000000000",
                model="mattersim_1_0_0_5m",
            )

    def test_list_and_filter_operations(self, test_project_id):
        """Test list operations and filtering"""

        # Test basic listing
        list_result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=10)
        
        assert "results" in list_result
        assert isinstance(list_result["results"], list)
        
        print(f"Listed {len(list_result['results'])} explorations")

        # Test pagination
        paginated_result = list_ht_mlrelax_explorations(
            project_id=test_project_id, limit=5, offset=0
        )
        
        assert "results" in paginated_result
        assert len(paginated_result["results"]) <= 5

        # If there are results, test structure
        if list_result["results"]:
            exploration = list_result["results"][0]
            expected_fields = ["id", "name", "model", "source_ht_sqs_exploration"]
            
            for field in expected_fields:
                assert field in exploration, f"Missing field: {field}"
                
            print(f"Sample exploration structure: {exploration}")
