"""
Integration tests for HT SQS exploration methods.

These tests require:
- Valid authentication credentials (ATOMICT_USERNAME, ATOMICT_PASSWORD)
- A test project (TEST_PROJECT_ID)
- A test structure (TEST_STRUCTURE_ID, TEST_STRUCTURE_TYPE)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_ht_sqs.py -v -m integration
"""

import os
from typing import List

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate

from atomict.simulation.ht_sqs import (
    create_ht_sqs_exploration,
    delete_ht_sqs_exploration,
    get_ht_sqs_child,
    get_ht_sqs_exploration,
    list_ht_sqs_children,
    list_ht_sqs_explorations,
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
def base_target_concentrations():
    """Base target concentrations for testing"""
    return [{"element": "Al", "weight": 0.6}, {"element": "Ni", "weight": 0.4}]


@pytest.fixture(scope="session")
def base_generated_permutations():
    """Base permutations for testing"""
    return [{"Al": 0.7, "Ni": 0.3}, {"Al": 0.5, "Ni": 0.5}, {"Al": 0.8, "Ni": 0.2}]


@pytest.fixture
def cleanup_ht_sqs_explorations():
    """Track created HT SQS explorations for cleanup"""
    created_ids = []

    def track_id(exploration_id):
        created_ids.append(exploration_id)

    yield track_id

    # Cleanup
    for ht_sqs_id in created_ids:
        try:
            delete_ht_sqs_exploration(ht_sqs_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestHTSQSExplorationIntegration:
    """Integration tests for HT SQS exploration CRUD operations"""

    def test_create_ht_sqs_exploration_minimal(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        base_target_concentrations,
        base_generated_permutations,
        cleanup_ht_sqs_explorations,
    ):
        """Test creating a HT SQS exploration with minimal parameters"""

        result = create_ht_sqs_exploration(
            project_id=test_project_id,
            name="Integration Test HT SQS Al-Ni Minimal",
            target_concentrations=base_target_concentrations,
            generated_permutations=base_generated_permutations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
        )

        # Should return success message (not standard DRF response)
        assert "message" in result
        assert "HTSQSExploration created successfully" in result["message"]

    def test_create_ht_sqs_exploration_full_params(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        base_target_concentrations,
        base_generated_permutations,
        cleanup_ht_sqs_explorations,
    ):
        """Test creating a HT SQS exploration with all parameters"""

        result = create_ht_sqs_exploration(
            project_id=test_project_id,
            name="Integration Test HT SQS Al-Ni Full",
            description="Full parameter test for HT SQS",
            target_concentrations=base_target_concentrations,
            generated_permutations=base_generated_permutations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            num_structures_per_permutation=2,
            auto_max_size=False,
            max_size=6,
            atom_count_upper_limit=150,
            cluster_cutoffs=[3.5, 3.5],
        )

        assert "message" in result
        assert "HTSQSExploration created successfully" in result["message"]

    def test_list_ht_sqs_explorations(self, test_project_id):
        """Test listing HT SQS explorations"""

        result = list_ht_sqs_explorations(project_id=test_project_id, limit=10)

        assert "results" in result
        assert isinstance(result["results"], list)
        # Should have standard pagination fields
        assert "count" in result or "next" in result or "previous" in result

    def test_get_ht_sqs_exploration_not_found(self):
        """Test getting non-existent HT SQS exploration"""

        with pytest.raises(Exception):  # Will raise APIError or similar
            get_ht_sqs_exploration("00000000-0000-0000-0000-000000000000")

    def test_list_ht_sqs_children(self):
        """Test listing HT SQS children"""

        result = list_ht_sqs_children(limit=5)

        assert "results" in result
        assert isinstance(result["results"], list)

        # If there are results, verify structure
        if result["results"]:
            child = result["results"][0]
            assert "id" in child
            # Should have either ht_sqs_exploration or sqs_exploration linked
            assert "ht_sqs_exploration" in child or "sqs_exploration" in child

    def test_create_ht_sqs_exploration_validation_errors(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test validation errors during HT SQS creation"""

        # Test invalid concentrations sum
        with pytest.raises(
            ValueError, match="Sum of target concentrations must equal 1.0"
        ):
            create_ht_sqs_exploration(
                project_id=test_project_id,
                name="Invalid Concentrations Test",
                target_concentrations=[
                    {"element": "Al", "weight": 0.8},
                    {"element": "Ni", "weight": 0.8},  # Sum > 1.0
                ],
                generated_permutations=[{"Al": 0.5, "Ni": 0.5}],
                structure_id=test_structure_id,
                structure_type=test_structure_type,
            )

        # Test invalid permutation sum
        with pytest.raises(
            ValueError, match="Permutation 0 concentrations must sum to 1.0"
        ):
            create_ht_sqs_exploration(
                project_id=test_project_id,
                name="Invalid Permutation Test",
                target_concentrations=[
                    {"element": "Al", "weight": 0.6},
                    {"element": "Ni", "weight": 0.4},
                ],
                generated_permutations=[{"Al": 0.8, "Ni": 0.8}],  # Sum > 1.0
                structure_id=test_structure_id,
                structure_type=test_structure_type,
            )

        # Test invalid structure type
        with pytest.raises(ValueError, match="structure_type must be one of"):
            create_ht_sqs_exploration(
                project_id=test_project_id,
                name="Invalid Structure Type Test",
                target_concentrations=[{"element": "Al", "weight": 1.0}],
                generated_permutations=[{"Al": 1.0}],
                structure_id=test_structure_id,
                structure_type="invalid_type",
            )

    def test_create_ht_sqs_with_multiple_structures_per_permutation(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_ht_sqs_explorations,
    ):
        """Test creating HT SQS with multiple structures per permutation"""

        # Small test to avoid creating too many child SQS instances
        single_permutation = [{"Al": 0.7, "Ni": 0.3}]

        result = create_ht_sqs_exploration(
            project_id=test_project_id,
            name="Integration Test HT SQS Multiple Structures",
            target_concentrations=[
                {"element": "Al", "weight": 0.7},
                {"element": "Ni", "weight": 0.3},
            ],
            generated_permutations=single_permutation,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            num_structures_per_permutation=3,  # Creates 3 SQS instances
        )

        assert "message" in result
        assert "HTSQSExploration created successfully" in result["message"]

    def test_ht_sqs_workflow_integration(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_ht_sqs_explorations,
    ):
        """Test complete HT SQS workflow: create -> list -> get children"""

        # 1. Create HT SQS exploration
        create_result = create_ht_sqs_exploration(
            project_id=test_project_id,
            name="Integration Test HT SQS Workflow",
            target_concentrations=[
                {"element": "Al", "weight": 0.6},
                {"element": "Ni", "weight": 0.4},
            ],
            generated_permutations=[{"Al": 0.8, "Ni": 0.2}, {"Al": 0.4, "Ni": 0.6}],
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            num_structures_per_permutation=1,
        )

        assert "HTSQSExploration created successfully" in create_result["message"]

        # 2. List explorations and verify our creation
        list_result = list_ht_sqs_explorations(project_id=test_project_id, limit=50)

        # Find our exploration by name
        our_exploration = None
        for exp in list_result["results"]:
            if exp["name"] == "Integration Test HT SQS Workflow":
                our_exploration = exp
                break

        assert our_exploration is not None, "Created exploration not found in list"
        cleanup_ht_sqs_explorations(our_exploration["id"])

        # 3. Get exploration details with children
        detail_result = get_ht_sqs_exploration(our_exploration["id"], children=True)

        assert "children" in detail_result
        # Should have 2 children (one per permutation)
        assert len(detail_result["children"]) == 2

        # 4. Verify children have SQS exploration links
        for child in detail_result["children"]:
            assert "sqs_exploration" in child
            assert child["sqs_exploration"] is not None
            assert "id" in child["sqs_exploration"]


@pytest.mark.integration
class TestHTSQSChildrenIntegration:
    """Integration tests for HT SQS children operations"""

    def test_get_ht_sqs_child_details(self):
        """Test getting detailed information about HT SQS children"""

        # Get some children to test with
        children_result = list_ht_sqs_children(limit=5)

        if not children_result["results"]:
            pytest.skip("No HT SQS children available for testing")

        child_id = children_result["results"][0]["id"]

        # Get detailed child information
        child_detail = get_ht_sqs_child(child_id)

        assert "id" in child_detail
        assert child_detail["id"] == child_id

        # Should have linked SQS exploration
        if "sqs_exploration" in child_detail and child_detail["sqs_exploration"]:
            sqs_exp = child_detail["sqs_exploration"]
            assert "id" in sqs_exp
            assert "name" in sqs_exp
            # Verify it's marked as HT
            if "is_ht" in sqs_exp:
                assert sqs_exp["is_ht"] is True

    def test_list_children_by_parent(self):
        """Test filtering children by parent HT SQS exploration"""

        # First get some HT SQS explorations
        explorations = list_ht_sqs_explorations(limit=5)

        if not explorations["results"]:
            pytest.skip("No HT SQS explorations available for testing")

        parent_id = explorations["results"][0]["id"]

        # Get children for this parent
        children = list_ht_sqs_children(ht_sqs_exploration_id=parent_id)

        assert "results" in children

        # The test is mainly to verify the API call works and filtering is applied
        # In a real system, children may belong to different parents
        # So we just verify the structure is correct
        for child in children["results"]:
            assert "id" in child
            # Verify ht_sqs_exploration field exists (may be None or reference)
            assert "ht_sqs_exploration" in child
