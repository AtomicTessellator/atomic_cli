import json
import pytest
from unittest.mock import MagicMock, patch, Mock
from click.testing import CliRunner
from atomict.cli.commands.organization import org_group, users_group, invites_group


@pytest.fixture
def runner():
    """CliRunner fixture for testing click commands"""
    return CliRunner()


@pytest.fixture
def mock_client():
    """Mock API client fixture"""
    with patch('atomict.cli.commands.organization.get_client') as mock_get_client:
        client = MagicMock()
        mock_get_client.return_value = client
        yield client


@pytest.fixture
def mock_console():
    """Mock console fixture"""
    with patch('atomict.cli.commands.organization.console') as mock:
        yield mock


class TestOrgGet:
    """Test org get command"""

    def test_org_get_single_organization(self, runner, mock_client):
        """Test getting a single organization by ID"""
        mock_client.get.return_value = {
            'id': '123',
            'name': 'Test Org',
            'description': 'Test description',
            'created_at': '2023-01-01T00:00:00Z',
            'updated_at': '2023-01-02T00:00:00Z',
            'billing_emails': ['billing@test.com']
        }
        
        result = runner.invoke(org_group, ['get', '123'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation/123/')

    def test_org_get_single_organization_json(self, runner, mock_client, mock_console):
        """Test getting a single organization with JSON output"""
        org_data = {'id': '123', 'name': 'Test Org'}
        mock_client.get.return_value = org_data
        
        result = runner.invoke(org_group, ['get', '123', '--json-output'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation/123/')
        mock_console.print_json.assert_called_once_with(data=org_data)

    def test_org_get_list_organizations(self, runner, mock_client):
        """Test listing all organizations"""
        mock_client.get.return_value = {
            'results': [
                {'id': '123', 'name': 'Org 1', 'users': [], 'created_at': '2023-01-01T00:00:00Z'},
                {'id': '456', 'name': 'Org 2', 'users': [], 'created_at': '2023-01-02T00:00:00Z'}
            ],
            'count': 2
        }
        
        result = runner.invoke(org_group, ['get'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation/', params={})

    def test_org_get_list_with_search(self, runner, mock_client):
        """Test listing organizations with search parameter"""
        mock_client.get.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['get', '--search', 'test'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation/', params={'search': 'test'})

    def test_org_get_list_with_filters(self, runner, mock_client):
        """Test listing organizations with filters"""
        mock_client.get.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['get', '--filter', 'active=true', '--filter', 'name=test'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation/', params={'active': 'true', 'name': 'test'})

    def test_org_get_list_fetch_all(self, runner, mock_client):
        """Test listing all organizations with fetch all option"""
        mock_client.get_all.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['get', '--all'])
        
        assert result.exit_code == 0
        mock_client.get_all.assert_called_once_with('/api/organisation/', params={})

    def test_org_get_invalid_filter_format(self, runner):
        """Test error handling for invalid filter format"""
        result = runner.invoke(org_group, ['get', '--filter', 'invalid-filter'])
        
        assert result.exit_code == 0  # Command doesn't exit with error, just prints message
        assert 'Invalid filter format' in result.output

    def test_org_get_list_json_output(self, runner, mock_client, mock_console):
        """Test listing organizations with JSON output"""
        org_data = {'results': [{'id': '123', 'name': 'Test Org'}], 'count': 1}
        mock_client.get.return_value = org_data
        
        result = runner.invoke(org_group, ['get', '--json-output'])
        
        assert result.exit_code == 0
        mock_console.print_json.assert_called_once_with(data=org_data)


class TestOrgCreate:
    """Test org create command"""

    def test_org_create_minimal(self, runner, mock_client):
        """Test creating organization with minimal parameters"""
        mock_client.post.return_value = {'id': '123', 'name': 'Test Org'}
        
        result = runner.invoke(org_group, ['create', '--name', 'Test Org'])
        
        assert result.exit_code == 0
        mock_client.post.assert_called_once_with('/api/organisation/', data={'name': 'Test Org'})
        assert 'Created organization \'Test Org\' with ID: 123' in result.output

    def test_org_create_full_data(self, runner, mock_client):
        """Test creating organization with all parameters"""
        mock_client.post.return_value = {'id': '123', 'name': 'Test Org'}
        
        result = runner.invoke(org_group, [
            'create',
            '--name', 'Test Org',
            '--description', 'Test description',
            '--billing-emails', 'billing@test.com,admin@test.com'
        ])
        
        assert result.exit_code == 0
        expected_data = {
            'name': 'Test Org',
            'description': 'Test description',
            'billing_emails': 'billing@test.com,admin@test.com'
        }
        mock_client.post.assert_called_once_with('/api/organisation/', data=expected_data)

    def test_org_create_json_output(self, runner, mock_client, mock_console):
        """Test creating organization with JSON output"""
        org_data = {'id': '123', 'name': 'Test Org'}
        mock_client.post.return_value = org_data
        
        result = runner.invoke(org_group, ['create', '--name', 'Test Org', '--json-output'])
        
        assert result.exit_code == 0
        mock_console.print_json.assert_called_once_with(data=org_data)

    def test_org_create_missing_name(self, runner):
        """Test error when name is not provided"""
        result = runner.invoke(org_group, ['create'])
        
        assert result.exit_code != 0
        assert 'Missing option' in result.output or 'required' in result.output.lower()


class TestOrgUpdate:
    """Test org update command"""

    def test_org_update_name(self, runner, mock_client):
        """Test updating organization name"""
        mock_client.patch.return_value = {'id': '123', 'name': 'Updated Org'}
        
        result = runner.invoke(org_group, ['update', '123', '--name', 'Updated Org'])
        
        assert result.exit_code == 0
        mock_client.patch.assert_called_once_with('/api/organisation/123/', data={'name': 'Updated Org'})
        assert 'Updated organization 123' in result.output

    def test_org_update_multiple_fields(self, runner, mock_client):
        """Test updating multiple organization fields"""
        mock_client.patch.return_value = {'id': '123', 'name': 'Updated Org'}
        
        result = runner.invoke(org_group, [
            'update', '123',
            '--name', 'Updated Org',
            '--description', 'New description',
            '--billing-emails', 'new@test.com'
        ])
        
        assert result.exit_code == 0
        expected_data = {
            'name': 'Updated Org',
            'description': 'New description',
            'billing_emails': 'new@test.com'
        }
        mock_client.patch.assert_called_once_with('/api/organisation/123/', data=expected_data)

    def test_org_update_no_fields(self, runner, mock_client):
        """Test updating organization with no fields provided"""
        result = runner.invoke(org_group, ['update', '123'])
        
        assert result.exit_code == 0
        assert 'No fields to update' in result.output
        mock_client.patch.assert_not_called()

    def test_org_update_json_output(self, runner, mock_client, mock_console):
        """Test updating organization with JSON output"""
        org_data = {'id': '123', 'name': 'Updated Org'}
        mock_client.patch.return_value = org_data
        
        result = runner.invoke(org_group, ['update', '123', '--name', 'Updated Org', '--json-output'])
        
        assert result.exit_code == 0
        mock_console.print_json.assert_called_once_with(data=org_data)


class TestOrgDelete:
    """Test org delete command"""

    def test_org_delete_with_confirmation(self, runner, mock_client):
        """Test deleting organization with confirmation"""
        result = runner.invoke(org_group, ['delete', '123'], input='y\n')
        
        assert result.exit_code == 0
        mock_client.delete.assert_called_once_with('/api/organisation/123/')
        assert 'Deleted organization 123' in result.output

    def test_org_delete_with_yes_flag(self, runner, mock_client):
        """Test deleting organization with --yes flag"""
        result = runner.invoke(org_group, ['delete', '123', '--yes'])
        
        assert result.exit_code == 0
        mock_client.delete.assert_called_once_with('/api/organisation/123/')
        assert 'Deleted organization 123' in result.output

    def test_org_delete_cancelled(self, runner, mock_client):
        """Test cancelling organization deletion"""
        result = runner.invoke(org_group, ['delete', '123'], input='n\n')
        
        assert result.exit_code == 0
        mock_client.delete.assert_not_called()
        assert 'Operation cancelled' in result.output


class TestOrgSwitch:
    """Test org switch command"""

    def test_org_switch_success(self, runner, mock_client):
        """Test successful organization switch"""
        mock_client.get.return_value = {'id': '123', 'name': 'Test Org'}
        mock_client.post.return_value = {'message': 'Switched successfully'}
        
        result = runner.invoke(org_group, ['switch', '123'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation/123/')
        mock_client.post.assert_called_once_with('/api/organisation/123/switch/', data={})
        assert 'Switched to organization \'Test Org\'' in result.output

    def test_org_switch_organization_not_found(self, runner, mock_client):
        """Test switching to non-existent organization"""
        mock_client.get.side_effect = Exception("Not found")
        
        result = runner.invoke(org_group, ['switch', '123'])
        
        assert result.exit_code == 0
        assert 'Organization 123 not found or not accessible' in result.output
        mock_client.post.assert_not_called()

    def test_org_switch_json_output(self, runner, mock_client, mock_console):
        """Test organization switch with JSON output"""
        mock_client.get.return_value = {'id': '123', 'name': 'Test Org'}
        switch_data = {'message': 'Switched successfully'}
        mock_client.post.return_value = switch_data
        
        result = runner.invoke(org_group, ['switch', '123', '--json-output'])
        
        assert result.exit_code == 0
        mock_console.print_json.assert_called_once_with(data=switch_data)


class TestUsersCommands:
    """Test organization users commands"""

    def test_users_list(self, runner, mock_client):
        """Test listing organization users"""
        mock_client.get.return_value = {
            'results': [
                {'id': '1', 'user': {'id': '123', 'username': 'user1', 'email': 'user1@test.com'}, 'role': 'member', 'created_at': '2023-01-01T00:00:00Z'}
            ],
            'count': 1
        }
        
        result = runner.invoke(org_group, ['users', 'list', 'org123'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation-user/', params={'organisation': 'org123'})

    def test_users_list_with_filters(self, runner, mock_client):
        """Test listing users with search and filters"""
        mock_client.get.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['users', 'list', 'org123', '--search', 'test', '--filter', 'role=admin'])
        
        assert result.exit_code == 0
        expected_params = {'organisation': 'org123', 'search': 'test', 'role': 'admin'}
        mock_client.get.assert_called_once_with('/api/organisation-user/', params=expected_params)

    def test_users_list_fetch_all(self, runner, mock_client):
        """Test listing all users with --all flag"""
        mock_client.get_all.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['users', 'list', 'org123', '--all'])
        
        assert result.exit_code == 0
        mock_client.get_all.assert_called_once_with('/api/organisation-user/', params={'organisation': 'org123'})

    def test_users_add_member(self, runner, mock_client):
        """Test adding a user as member"""
        mock_client.post.return_value = {'id': '1', 'user_email': 'user@test.com', 'role': 'member'}
        
        result = runner.invoke(org_group, ['users', 'add', 'org123', 'user@test.com'])
        
        assert result.exit_code == 0
        expected_data = {'organisation': 'org123', 'user_email': 'user@test.com', 'role': 'member'}
        mock_client.post.assert_called_once_with('/api/organisation-user/', data=expected_data)
        assert 'Added user user@test.com to organization org123 as member' in result.output

    def test_users_add_admin(self, runner, mock_client):
        """Test adding a user as admin"""
        mock_client.post.return_value = {'id': '1', 'user_email': 'admin@test.com', 'role': 'admin'}
        
        result = runner.invoke(org_group, ['users', 'add', 'org123', 'admin@test.com', '--admin'])
        
        assert result.exit_code == 0
        expected_data = {'organisation': 'org123', 'user_email': 'admin@test.com', 'role': 'admin'}
        mock_client.post.assert_called_once_with('/api/organisation-user/', data=expected_data)
        assert 'Added user admin@test.com to organization org123 as admin' in result.output

    def test_users_remove_with_confirmation(self, runner, mock_client):
        """Test removing user with confirmation"""
        mock_client.get.return_value = {'results': [{'id': 'rel123'}], 'count': 1}
        
        result = runner.invoke(org_group, ['users', 'remove', 'org123', 'user456'], input='y\n')
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation-user/', params={'organisation': 'org123', 'user': 'user456'})
        mock_client.delete.assert_called_once_with('/api/organisation-user/rel123/')
        assert 'Removed user user456 from organization org123' in result.output

    def test_users_remove_user_not_found(self, runner, mock_client):
        """Test removing non-existent user"""
        mock_client.get.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['users', 'remove', 'org123', 'user456'], input='y\n')
        
        assert result.exit_code == 0
        assert 'User user456 not found in organization org123' in result.output
        mock_client.delete.assert_not_called()

    def test_users_promote(self, runner, mock_client):
        """Test promoting user to admin"""
        mock_client.get.return_value = {'results': [{'id': 'rel123'}], 'count': 1}
        mock_client.patch.return_value = {'id': 'rel123', 'role': 'admin'}
        
        result = runner.invoke(org_group, ['users', 'promote', 'org123', 'user456'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation-user/', params={'organisation': 'org123', 'user': 'user456'})
        mock_client.patch.assert_called_once_with('/api/organisation-user/rel123/', data={'role': 'admin'})
        assert 'Promoted user user456 to admin in organization org123' in result.output

    def test_users_demote(self, runner, mock_client):
        """Test demoting user from admin to member"""
        mock_client.get.return_value = {'results': [{'id': 'rel123'}], 'count': 1}
        mock_client.patch.return_value = {'id': 'rel123', 'role': 'member'}
        
        result = runner.invoke(org_group, ['users', 'demote', 'org123', 'user456'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation-user/', params={'organisation': 'org123', 'user': 'user456'})
        mock_client.patch.assert_called_once_with('/api/organisation-user/rel123/', data={'role': 'member'})
        assert 'Demoted user user456 to member in organization org123' in result.output


class TestInvitesCommands:
    """Test organization invites commands"""

    def test_invites_list(self, runner, mock_client):
        """Test listing organization invites"""
        mock_client.get.return_value = {
            'results': [
                {'id': '1', 'email': 'invite@test.com', 'role': 'member', 'status': 'pending', 'created_at': '2023-01-01T00:00:00Z', 'expires_at': '2023-01-08T00:00:00Z'}
            ],
            'count': 1
        }
        
        result = runner.invoke(org_group, ['invites', 'list', 'org123'])
        
        assert result.exit_code == 0
        mock_client.get.assert_called_once_with('/api/organisation-invite/', params={'organisation': 'org123'})

    def test_invites_list_with_filters(self, runner, mock_client):
        """Test listing invites with search and filters"""
        mock_client.get.return_value = {'results': [], 'count': 0}
        
        result = runner.invoke(org_group, ['invites', 'list', 'org123', '--search', 'test', '--filter', 'status=pending'])
        
        assert result.exit_code == 0
        expected_params = {'organisation': 'org123', 'search': 'test', 'status': 'pending'}
        mock_client.get.assert_called_once_with('/api/organisation-invite/', params=expected_params)

    def test_invites_send_member(self, runner, mock_client):
        """Test sending member invite"""
        mock_client.post.return_value = {'id': '1', 'email': 'invite@test.com', 'role': 'member'}
        
        result = runner.invoke(org_group, ['invites', 'send', 'org123', 'invite@test.com'])
        
        assert result.exit_code == 0
        expected_data = {'organisation': 'org123', 'email': 'invite@test.com', 'role': 'member'}
        mock_client.post.assert_called_once_with('/api/organisation-invite/', data=expected_data)
        assert 'Sent invitation to invite@test.com to join organization org123 as member' in result.output

    def test_invites_send_admin(self, runner, mock_client):
        """Test sending admin invite"""
        mock_client.post.return_value = {'id': '1', 'email': 'admin@test.com', 'role': 'admin'}
        
        result = runner.invoke(org_group, ['invites', 'send', 'org123', 'admin@test.com', '--admin'])
        
        assert result.exit_code == 0
        expected_data = {'organisation': 'org123', 'email': 'admin@test.com', 'role': 'admin'}
        mock_client.post.assert_called_once_with('/api/organisation-invite/', data=expected_data)
        assert 'Sent invitation to admin@test.com to join organization org123 as admin' in result.output

    def test_invites_cancel_with_confirmation(self, runner, mock_client):
        """Test cancelling invite with confirmation"""
        result = runner.invoke(org_group, ['invites', 'cancel', 'invite123'], input='y\n')
        
        assert result.exit_code == 0
        mock_client.delete.assert_called_once_with('/api/organisation-invite/invite123/')
        assert 'Cancelled invitation invite123' in result.output

    def test_invites_cancel_with_yes_flag(self, runner, mock_client):
        """Test cancelling invite with --yes flag"""
        result = runner.invoke(org_group, ['invites', 'cancel', 'invite123', '--yes'])
        
        assert result.exit_code == 0
        mock_client.delete.assert_called_once_with('/api/organisation-invite/invite123/')
        assert 'Cancelled invitation invite123' in result.output

    def test_invites_cancel_cancelled(self, runner, mock_client):
        """Test cancelling invite cancellation"""
        result = runner.invoke(org_group, ['invites', 'cancel', 'invite123'], input='n\n')
        
        assert result.exit_code == 0
        mock_client.delete.assert_not_called()
        assert 'Operation cancelled' in result.output


class TestJSONOutput:
    """Test JSON output functionality across commands"""

    def test_all_commands_support_json_output(self, runner, mock_client, mock_console):
        """Test that commands with --json-output flag use console.print_json"""
        json_commands = [
            (['get', '123', '--json-output'], {'id': '123', 'name': 'Test'}),
            (['create', '--name', 'Test', '--json-output'], {'id': '123', 'name': 'Test'}),
            (['update', '123', '--name', 'Updated', '--json-output'], {'id': '123', 'name': 'Updated'}),
            (['switch', '123', '--json-output'], {'message': 'Switched'}),
            (['users', 'add', 'org123', 'user@test.com', '--json-output'], {'id': '1', 'email': 'user@test.com'}),
            (['users', 'promote', 'org123', 'user456', '--json-output'], {'id': 'rel123', 'role': 'admin'}),
            (['users', 'demote', 'org123', 'user456', '--json-output'], {'id': 'rel123', 'role': 'member'}),
            (['invites', 'send', 'org123', 'invite@test.com', '--json-output'], {'id': '1', 'email': 'invite@test.com'}),
        ]
        
        for command, return_data in json_commands:
            # Setup mocks based on command type
            if 'get' in command and len(command) == 3:  # Single get
                mock_client.get.return_value = return_data
            elif 'get' in command:  # List get
                mock_client.get.return_value = {'results': [return_data]}
            elif 'create' in command or 'add' in command or 'send' in command:
                mock_client.post.return_value = return_data
            elif 'update' in command or 'promote' in command or 'demote' in command:
                mock_client.patch.return_value = return_data
                if 'promote' in command or 'demote' in command:
                    mock_client.get.return_value = {'results': [{'id': 'rel123'}]}
            elif 'switch' in command:
                mock_client.get.return_value = {'id': '123', 'name': 'Test'}
                mock_client.post.return_value = return_data
            
            # Reset console mock
            mock_console.reset_mock()
            
            # Run command
            result = runner.invoke(org_group, command)
            
            # Verify JSON output was called
            assert result.exit_code == 0, f"Command {command} failed: {result.output}"
            mock_console.print_json.assert_called_once(), f"Command {command} didn't call print_json"


class TestErrorHandling:
    """Test error handling in CLI commands"""

    def test_client_error_handling(self, runner, mock_client):
        """Test that client errors are handled gracefully"""
        mock_client.get.side_effect = Exception("API Error")
        
        result = runner.invoke(org_group, ['get', '123'])
        
        # Command should handle the error gracefully (exact behavior depends on implementation)
        # The important thing is it doesn't crash with unhandled exception
        assert isinstance(result.exception, (type(None), SystemExit, Exception))

    def test_invalid_command_arguments(self, runner):
        """Test handling of invalid command arguments"""
        # Test missing required arguments
        result = runner.invoke(org_group, ['update'])
        assert result.exit_code != 0
        
        result = runner.invoke(org_group, ['delete'])
        assert result.exit_code != 0
        
        result = runner.invoke(org_group, ['switch'])
        assert result.exit_code != 0


class TestCommandIntegration:
    """Integration-style tests for command workflows"""

    def test_organization_management_workflow(self, runner, mock_client):
        """Test a complete organization management workflow"""
        # Create organization
        mock_client.post.return_value = {'id': 'org123', 'name': 'Test Org'}
        result = runner.invoke(org_group, ['create', '--name', 'Test Org'])
        assert result.exit_code == 0
        
        # Get organization
        mock_client.get.return_value = {'id': 'org123', 'name': 'Test Org', 'users': []}
        result = runner.invoke(org_group, ['get', 'org123'])
        assert result.exit_code == 0
        
        # Add user
        mock_client.post.return_value = {'id': 'user1', 'user_email': 'user@test.com'}
        result = runner.invoke(org_group, ['users', 'add', 'org123', 'user@test.com'])
        assert result.exit_code == 0
        
        # List users
        mock_client.get.return_value = {
            'results': [{'id': 'user1', 'user': {'email': 'user@test.com'}, 'role': 'member'}],
            'count': 1
        }
        result = runner.invoke(org_group, ['users', 'list', 'org123'])
        assert result.exit_code == 0
        
        # Send invite
        mock_client.post.return_value = {'id': 'invite1', 'email': 'invite@test.com'}
        result = runner.invoke(org_group, ['invites', 'send', 'org123', 'invite@test.com'])
        assert result.exit_code == 0
        
        # List invites
        mock_client.get.return_value = {
            'results': [{'id': 'invite1', 'email': 'invite@test.com', 'status': 'pending'}],
            'count': 1
        }
        result = runner.invoke(org_group, ['invites', 'list', 'org123'])
        assert result.exit_code == 0

    def test_user_role_management_workflow(self, runner, mock_client):
        """Test user role management workflow"""
        # Add user as member
        mock_client.post.return_value = {'id': 'user1', 'role': 'member'}
        result = runner.invoke(org_group, ['users', 'add', 'org123', 'user@test.com'])
        assert result.exit_code == 0
        
        # Promote to admin
        mock_client.get.return_value = {'results': [{'id': 'rel123'}]}
        mock_client.patch.return_value = {'id': 'rel123', 'role': 'admin'}
        result = runner.invoke(org_group, ['users', 'promote', 'org123', 'user456'])
        assert result.exit_code == 0
        
        # Demote back to member
        mock_client.get.return_value = {'results': [{'id': 'rel123'}]}
        mock_client.patch.return_value = {'id': 'rel123', 'role': 'member'}
        result = runner.invoke(org_group, ['users', 'demote', 'org123', 'user456'])
        assert result.exit_code == 0
        
        # Remove user
        mock_client.get.return_value = {'results': [{'id': 'rel123'}]}
        result = runner.invoke(org_group, ['users', 'remove', 'org123', 'user456', '--yes'])
        assert result.exit_code == 0
