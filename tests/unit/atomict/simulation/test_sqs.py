import pytest
from unittest.mock import Mock, patch
from atomict.simulation.sqs import (
    get_simulation,
    create_sqs_exploration,
    delete_sqs_exploration,
    associate_user_upload_with_sqs_simulation,
)


@pytest.fixture
def mock_get():
    with patch('atomict.simulation.sqs.get') as mock:
        yield mock


@pytest.fixture
def mock_post():
    with patch('atomict.simulation.sqs.post') as mock:
        yield mock


@pytest.fixture
def mock_delete():
    with patch('atomict.simulation.sqs.delete') as mock:
        yield mock


class TestGetSimulation:
    def test_get_simulation_without_params(self, mock_get):
        """Test getting a SQS simulation without additional parameters"""
        mock_get.return_value = {'id': '123', 'name': 'test_sqs'}
        
        result = get_simulation('123')
        
        assert result == {'id': '123', 'name': 'test_sqs'}
        mock_get.assert_called_once_with('api/sqs-exploration/123/')

    def test_get_simulation_with_full_param(self, mock_get):
        """Test getting a SQS simulation with full parameter"""
        mock_get.return_value = {'id': '123', 'name': 'test_sqs', 'details': 'full'}
        
        result = get_simulation('123', full=True)
        
        assert result == {'id': '123', 'name': 'test_sqs', 'details': 'full'}
        mock_get.assert_called_once_with('api/sqs-exploration/123/?full=true')

    def test_get_simulation_with_params(self, mock_get):
        """Test getting a SQS simulation with additional parameters"""
        mock_get.return_value = {'id': '123', 'name': 'test_sqs'}
        
        result = get_simulation('123', full=True, status='completed')
        
        assert result == {'id': '123', 'name': 'test_sqs'}
        # Query string parameter order can vary, so just check the call was made
        mock_get.assert_called_once()
        call_url = mock_get.call_args.args[0]
        assert call_url.startswith('api/sqs-exploration/123/?')
        assert 'full=true' in call_url
        assert 'status=completed' in call_url


class TestCreateSQSExploration:
    def test_create_sqs_exploration_minimal_params(self, mock_post):
        """Test creating a SQS exploration with minimal required parameters (defaults to userupload)"""
        mock_post.return_value = {'id': 'new_sqs_123'}
        target_concentrations = [
            {"element": "Cu", "weight": 0.5},
            {"element": "Zn", "weight": 0.5}
        ]
        
        result = create_sqs_exploration(
            project_id='proj_456',
            name='Test SQS',
            target_concentrations=target_concentrations,
            structure_id='upload_789'
        )
        
        assert result == {'id': 'new_sqs_123'}
        mock_post.assert_called_once_with(
            'api/sqs-exploration/',
            {
                'project': 'proj_456',
                'name': 'Test SQS',
                'description': None,
                'target_concentrations': target_concentrations,
                'action': 'DRAFT',
                'max_size': 8,
                'cluster_cutoffs': [4.0, 4.0],
                'starting_structure_userupload_id': 'upload_789'
            },
            extra_headers={'Content-Type': 'application/json'}
        )

    def test_create_sqs_exploration_fhiaims_structure(self, mock_post):
        """Test creating SQS exploration with FHI-aims starting structure"""
        mock_post.return_value = {'id': 'new_sqs_456'}
        target_concentrations = [
            {"element": "Fe", "weight": 0.25},
            {"element": "Ni", "weight": 0.25},
            {"element": "Cr", "weight": 0.5}
        ]
        extra_kwargs = {'custom_param': 'value'}
        
        result = create_sqs_exploration(
            project_id='proj_123',
            name='FHI-aims SQS Test',
            target_concentrations=target_concentrations,
            structure_id='fhiaims_456',
            structure_type='fhiaims',
            action='LAUNCH',
            description='Test description',
            max_size=16,
            cluster_cutoffs=[3.5, 4.5],
            extra_exploration_kwargs=extra_kwargs
        )
        
        assert result == {'id': 'new_sqs_456'}
        mock_post.assert_called_once_with(
            'api/sqs-exploration/',
            {
                'project': 'proj_123',
                'name': 'FHI-aims SQS Test',
                'description': 'Test description',
                'target_concentrations': target_concentrations,
                'action': 'LAUNCH',
                'max_size': 16,
                'cluster_cutoffs': [3.5, 4.5],
                'starting_structure_id': 'fhiaims_456',
                'custom_param': 'value'
            },
            extra_headers={'Content-Type': 'application/json'}
        )

    def test_create_sqs_exploration_mlrelax_structure(self, mock_post):
        """Test creating SQS exploration with MLRelax starting structure"""
        mock_post.return_value = {'id': 'new_sqs_789'}
        target_concentrations = [
            {"element": "Al", "weight": 0.7},
            {"element": "Mg", "weight": 0.3}
        ]
        
        result = create_sqs_exploration(
            project_id='proj_789',
            name='MLRelax SQS',
            target_concentrations=target_concentrations,
            structure_id='mlrelax_123',
            structure_type='mlrelax'
        )
        
        assert result == {'id': 'new_sqs_789'}
        mock_post.assert_called_once_with(
            'api/sqs-exploration/',
            {
                'project': 'proj_789',
                'name': 'MLRelax SQS',
                'description': None,
                'target_concentrations': target_concentrations,
                'action': 'DRAFT',
                'max_size': 8,
                'cluster_cutoffs': [4.0, 4.0],
                'starting_structure_mlrelax_id': 'mlrelax_123'
            },
            extra_headers={'Content-Type': 'application/json'}
        )

    def test_create_sqs_exploration_valid_structure_types(self, mock_post):
        """Test creating SQS explorations with all valid structure types"""
        target_concentrations = [{"element": "Cu", "weight": 1.0}]
        
        structure_types = [
            ('fhiaims', 'starting_structure_id'),
            ('mlrelax', 'starting_structure_mlrelax_id'),
            ('userupload', 'starting_structure_userupload_id')
        ]
        
        for structure_type, expected_field in structure_types:
            mock_post.return_value = {'id': f'sqs_{structure_type}'}
            
            result = create_sqs_exploration(
                project_id='proj_123',
                name=f'Test {structure_type}',
                target_concentrations=target_concentrations,
                structure_id='struct_456',
                structure_type=structure_type
            )
            
            assert result == {'id': f'sqs_{structure_type}'}
            # Check that the correct field was used
            call_payload = mock_post.call_args[0][1]
            assert expected_field in call_payload
            assert call_payload[expected_field] == 'struct_456'
            
            mock_post.reset_mock()

    def test_create_sqs_exploration_invalid_action(self, mock_post):
        """Test error handling for invalid action"""
        target_concentrations = [{"element": "Cu", "weight": 1.0}]
        
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_sqs_exploration(
                project_id='proj_123',
                name='Invalid Action',
                target_concentrations=target_concentrations,
                structure_id='struct_456',
                action='INVALID'
            )
        
        mock_post.assert_not_called()

    def test_create_sqs_exploration_invalid_structure_type(self, mock_post):
        """Test error handling for invalid structure type"""
        target_concentrations = [{"element": "Cu", "weight": 1.0}]
        
        with pytest.raises(ValueError, match="structure_type must be one of"):
            create_sqs_exploration(
                project_id='proj_123',
                name='Invalid Structure Type',
                target_concentrations=target_concentrations,
                structure_id='struct_456',
                structure_type='invalid_type'
            )
        
        mock_post.assert_not_called()

    def test_create_sqs_exploration_empty_target_concentrations(self, mock_post):
        """Test error handling for empty target concentrations"""
        with pytest.raises(ValueError, match="target_concentrations cannot be empty"):
            create_sqs_exploration(
                project_id='proj_123',
                name='Empty Concentrations',
                target_concentrations=[],
                structure_id='struct_456'
            )
        
        mock_post.assert_not_called()

    def test_create_sqs_exploration_invalid_concentration_format(self, mock_post):
        """Test error handling for invalid concentration format"""
        invalid_concentrations = [{"element": "Cu"}]  # Missing weight
        
        with pytest.raises(ValueError, match="Each target concentration must be a dict with 'element' and 'weight' keys"):
            create_sqs_exploration(
                project_id='proj_123',
                name='Invalid Format',
                target_concentrations=invalid_concentrations,
                structure_id='struct_456'
            )
        
        mock_post.assert_not_called()

    def test_create_sqs_exploration_invalid_weight_range(self, mock_post):
        """Test error handling for weight outside valid range"""
        invalid_concentrations = [{"element": "Cu", "weight": 1.5}]  # Weight > 1.0
        
        with pytest.raises(ValueError, match="Target concentration for Cu must be between 0 and 1"):
            create_sqs_exploration(
                project_id='proj_123',
                name='Invalid Weight',
                target_concentrations=invalid_concentrations,
                structure_id='struct_456'
            )
        
        mock_post.assert_not_called()

    def test_create_sqs_exploration_concentrations_not_sum_to_one(self, mock_post):
        """Test error handling when concentrations don't sum to 1.0"""
        invalid_concentrations = [
            {"element": "Cu", "weight": 0.3},
            {"element": "Zn", "weight": 0.3}  # Sum = 0.6, not 1.0
        ]
        
        with pytest.raises(ValueError, match="Sum of target concentrations must equal 1.0"):
            create_sqs_exploration(
                project_id='proj_123',
                name='Invalid Sum',
                target_concentrations=invalid_concentrations,
                structure_id='struct_456'
            )
        
        mock_post.assert_not_called()

    def test_create_sqs_exploration_concentrations_sum_tolerance(self, mock_post):
        """Test that small floating point errors in concentration sum are tolerated"""
        # Sum is 1.0000001, should be within tolerance
        concentrations = [
            {"element": "Cu", "weight": 0.3333333},
            {"element": "Zn", "weight": 0.3333333},
            {"element": "Al", "weight": 0.3333334}  # Sum = 1.0000000
        ]
        mock_post.return_value = {'id': 'tolerance_test'}
        
        result = create_sqs_exploration(
            project_id='proj_123',
            name='Tolerance Test',
            target_concentrations=concentrations,
            structure_id='struct_456'
        )
        
        assert result == {'id': 'tolerance_test'}
        mock_post.assert_called_once()


class TestDeleteSQSExploration:
    def test_delete_sqs_exploration(self, mock_delete):
        """Test deleting a SQS exploration"""
        mock_delete.return_value = {'status': 'deleted'}
        
        result = delete_sqs_exploration('sqs_123')
        
        assert result == {'status': 'deleted'}
        mock_delete.assert_called_once_with('api/sqs-exploration/sqs_123/')


class TestAssociateUserUploadWithSQSSimulation:
    def test_associate_user_upload_success(self, mock_post):
        """Test successful association of user upload with SQS simulation"""
        mock_post.return_value = {'id': 'association_123'}
        
        result = associate_user_upload_with_sqs_simulation('upload_456', 'sqs_789')
        
        assert result == {'id': 'association_123'}
        mock_post.assert_called_once_with(
            'api/sqs-simulation-file/',
            payload={'user_upload_id': 'upload_456', 'exploration_id': 'sqs_789'}
        )


class TestSQSModule:
    """Integration-style tests for the SQS module"""
    
    @patch('atomict.simulation.sqs.get')
    @patch('atomict.simulation.sqs.post')
    @patch('atomict.simulation.sqs.delete')
    def test_sqs_workflow_integration(self, mock_delete, mock_post, mock_get):
        """Test a typical SQS workflow using multiple functions"""
        # Mock responses for workflow
        target_concentrations = [
            {"element": "Cu", "weight": 0.5},
            {"element": "Zn", "weight": 0.5}
        ]
        
        mock_post.return_value = {'id': 'sqs_123'}
        mock_get.return_value = {'id': 'sqs_123', 'status': 'completed'}
        mock_delete.return_value = {'status': 'deleted'}
        
        # Create a SQS exploration
        create_result = create_sqs_exploration(
            project_id='proj_123',
            name='Integration Test SQS',
            target_concentrations=target_concentrations,
            structure_id='struct_456',
            structure_type='fhiaims',
            action='LAUNCH'
        )
        assert create_result['id'] == 'sqs_123'
        
        # Get the created exploration
        get_result = get_simulation('sqs_123')
        assert get_result['status'] == 'completed'
        
        # Delete the exploration
        delete_result = delete_sqs_exploration('sqs_123')
        assert delete_result['status'] == 'deleted'
        
        # Verify all API calls were made correctly
        mock_post.assert_called_once()
        mock_get.assert_called_once_with('api/sqs-exploration/sqs_123/')
        mock_delete.assert_called_once_with('api/sqs-exploration/sqs_123/')
