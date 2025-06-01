import pytest
from unittest.mock import Mock, patch
from atomict.simulation.phonon import (
    get_phonon_run,
    get_phono3py_sim_run,
    associate_user_upload_with_phono3py_sim_run,
    create_phono3py_run,
    get_phono3py_sim_run_files,
)
from atomict.simulation.models import (
    MODEL_ORB_D3_V2,
    MODEL_MATTERSIM_1_0_0_5M,
    MODEL_ORB_V3_CONSERVATIVE,
    MODEL_ESEN_30M_OAM,
)


@pytest.fixture
def mock_get():
    with patch('atomict.simulation.phonon.get') as mock:
        yield mock


@pytest.fixture
def mock_post():
    with patch('atomict.simulation.phonon.post') as mock:
        yield mock


class TestGetPhononRun:
    def test_get_phonon_run_without_params(self, mock_get):
        """Test getting a phonon run without additional parameters"""
        mock_get.return_value = {'id': '123', 'name': 'test_run'}
        
        result = get_phonon_run('123')
        
        assert result == {'id': '123', 'name': 'test_run'}
        mock_get.assert_called_once_with('api/phono3py-run/123/')

    def test_get_phonon_run_with_params(self, mock_get):
        """Test getting a phonon run with additional parameters"""
        mock_get.return_value = {'id': '123', 'name': 'test_run'}
        
        result = get_phonon_run('123', limit=10, offset=20)
        
        assert result == {'id': '123', 'name': 'test_run'}
        mock_get.assert_called_once_with('api/phono3py-run/123/?limit=10&offset=20')

    def test_get_phonon_run_with_single_param(self, mock_get):
        """Test getting a phonon run with a single parameter"""
        mock_get.return_value = {'id': '456'}
        
        result = get_phonon_run('456', status='completed')
        
        assert result == {'id': '456'}
        mock_get.assert_called_once_with('api/phono3py-run/456/?status=completed')


class TestGetPhono3pySimRun:
    def test_get_phono3py_sim_run_without_params(self, mock_get):
        """Test getting a phonon simulation run without additional parameters"""
        mock_get.return_value = {'id': '789', 'status': 'running'}
        
        result = get_phono3py_sim_run('789')
        
        assert result == {'id': '789', 'status': 'running'}
        mock_get.assert_called_once_with('api/phono3py-run-simulation/789/')

    def test_get_phono3py_sim_run_with_params(self, mock_get):
        """Test getting a phonon simulation run with additional parameters"""
        mock_get.return_value = {'id': '789', 'status': 'running'}
        
        result = get_phono3py_sim_run('789', include_files=True, expand_data=True)
        
        assert result == {'id': '789', 'status': 'running'}
        mock_get.assert_called_once_with('api/phono3py-run-simulation/789/?include_files=True&expand_data=True')


class TestAssociateUserUploadWithPhono3pySimRun:
    def test_associate_user_upload_success(self, mock_post):
        """Test successful association of user upload with phonon simulation run"""
        mock_post.return_value = {'id': 'association_123'}
        
        result = associate_user_upload_with_phono3py_sim_run('upload_456', 'run_789')
        
        assert result == {'id': 'association_123'}
        mock_post.assert_called_once_with(
            'api/phono3py-run-simulation-file/',
            payload={'user_upload_id': 'upload_456', 'phono3py_run_simulation_id': 'run_789'}
        )


class TestCreatePhono3pyRun:
    def test_create_phono3py_run_minimal_params(self, mock_post):
        """Test creating a phonon run with minimal required parameters"""
        mock_post.return_value = {'id': 'new_run_123'}
        
        result = create_phono3py_run(
            project_id='proj_456',
            source_geometry_id='geom_789',
            action='LAUNCH'
        )
        
        assert result == {'id': 'new_run_123'}
        mock_post.assert_called_once_with(
            'api/phono3py-run/',
            payload={
                'project_id': 'proj_456',
                'source_geometry_id': 'geom_789',
                'action': 'LAUNCH',
                'name': None,
                'description': None,
                'model': MODEL_ORB_D3_V2,
                'extra_simulation_kwargs': None,
            }
        )

    def test_create_phono3py_run_all_params(self, mock_post):
        """Test creating a phonon run with all parameters"""
        mock_post.return_value = {'id': 'new_run_456'}
        extra_kwargs = {'supercell': [2, 2, 2]}
        
        result = create_phono3py_run(
            project_id='proj_123',
            source_geometry_id='geom_456',
            action='DRAFT',
            name='Test Phonon Run',
            description='Test description',
            model=MODEL_MATTERSIM_1_0_0_5M,
            extra_simulation_kwargs=extra_kwargs
        )
        
        assert result == {'id': 'new_run_456'}
        mock_post.assert_called_once_with(
            'api/phono3py-run/',
            payload={
                'project_id': 'proj_123',
                'source_geometry_id': 'geom_456',
                'action': 'DRAFT',
                'name': 'Test Phonon Run',
                'description': 'Test description',
                'model': MODEL_MATTERSIM_1_0_0_5M,
                'extra_simulation_kwargs': extra_kwargs,
            }
        )

    def test_create_phono3py_run_valid_models(self, mock_post):
        """Test creating phonon runs with all valid model types"""
        valid_models = [
            MODEL_ORB_D3_V2,
            MODEL_MATTERSIM_1_0_0_5M,
            MODEL_ORB_V3_CONSERVATIVE,
            MODEL_ESEN_30M_OAM
        ]
        
        for model in valid_models:
            mock_post.return_value = {'id': f'run_{model}'}
            
            result = create_phono3py_run(
                project_id='proj_123',
                source_geometry_id='geom_456',
                action='LAUNCH',
                model=model
            )
            
            assert result == {'id': f'run_{model}'}

    def test_create_phono3py_run_invalid_model(self, mock_post):
        """Test error handling for invalid model"""
        with pytest.raises(ValueError, match="Invalid model: 999"):
            create_phono3py_run(
                project_id='proj_123',
                source_geometry_id='geom_456',
                action='LAUNCH',
                model=999
            )
        
        mock_post.assert_not_called()

    def test_create_phono3py_run_invalid_action(self, mock_post):
        """Test error handling for invalid action"""
        with pytest.raises(ValueError, match="Invalid action: INVALID"):
            create_phono3py_run(
                project_id='proj_123',
                source_geometry_id='geom_456',
                action='INVALID'
            )
        
        mock_post.assert_not_called()

    def test_create_phono3py_run_valid_actions(self, mock_post):
        """Test creating phonon runs with valid actions"""
        valid_actions = ['LAUNCH', 'DRAFT']
        
        for action in valid_actions:
            mock_post.return_value = {'id': f'run_{action.lower()}'}
            
            result = create_phono3py_run(
                project_id='proj_123',
                source_geometry_id='geom_456',
                action=action
            )
            
            assert result == {'id': f'run_{action.lower()}'}


class TestGetPhono3pySimRunFiles:
    def test_get_phono3py_sim_run_files(self, mock_get):
        """Test getting files associated with a phonon simulation run"""
        expected_files = [
            {'id': 'file_1', 'filename': 'POSCAR'},
            {'id': 'file_2', 'filename': 'OUTCAR'}
        ]
        mock_get.return_value = expected_files
        
        result = get_phono3py_sim_run_files('sim_run_123')
        
        assert result == expected_files
        mock_get.assert_called_once_with('api/phono3py-run-simulation-file/?phono3py_run_simulation__id=sim_run_123')

    def test_get_phono3py_sim_run_files_empty(self, mock_get):
        """Test getting files when no files are associated"""
        mock_get.return_value = []
        
        result = get_phono3py_sim_run_files('sim_run_456')
        
        assert result == []
        mock_get.assert_called_once_with('api/phono3py-run-simulation-file/?phono3py_run_simulation__id=sim_run_456')


class TestPhononModule:
    """Integration-style tests for the phonon module"""
    
    @patch('atomict.simulation.phonon.get')
    @patch('atomict.simulation.phonon.post')
    def test_phonon_workflow_integration(self, mock_post, mock_get):
        """Test a typical workflow using multiple phonon functions"""
        # Mock responses for workflow
        mock_post.return_value = {'id': 'run_123'}
        mock_get.side_effect = [
            {'id': 'run_123', 'status': 'completed'},  # get_phonon_run
            {'id': 'sim_123', 'status': 'completed'},  # get_phono3py_sim_run  
            [{'id': 'file_1', 'filename': 'phonopy.yaml'}]  # get_phono3py_sim_run_files
        ]
        
        # Create a phonon run
        run_result = create_phono3py_run(
            project_id='proj_123',
            source_geometry_id='geom_456',
            action='LAUNCH',
            name='Integration Test Run'
        )
        assert run_result['id'] == 'run_123'
        
        # Get the created run
        run_data = get_phonon_run('run_123')
        assert run_data['status'] == 'completed'
        
        # Get simulation data
        sim_data = get_phono3py_sim_run('sim_123')
        assert sim_data['status'] == 'completed'
        
        # Get associated files
        files = get_phono3py_sim_run_files('sim_123')
        assert len(files) == 1
        assert files[0]['filename'] == 'phonopy.yaml'
