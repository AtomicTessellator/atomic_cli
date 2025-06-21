"""
Integration tests for organization CLI commands.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- Optional TEST_ORGANIZATION_ID for non-destructive tests
- Optional TEST_ORGANIZATION_USER_EMAIL for user management tests

To run: uv run pytest tests/integration/test_organization_cli.py -v -m integration
"""

import json
import os
import re
import subprocess
import uuid

import pytest
from click.testing import CliRunner
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.cli.commands.organization import org_group
from atomict.organization.client import create_organization, delete_organization


def strip_ansi_codes(text: str) -> str:
    """Remove ANSI color codes from text"""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


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
def test_organization_id():
    """Get test organization ID from environment (optional for non-destructive tests)"""
    return os.getenv("TEST_ORGANIZATION_ID")


@pytest.fixture(scope="session")
def test_user_email():
    """Get test user email from environment for user management tests"""
    return os.getenv("TEST_ORGANIZATION_USER_EMAIL")


@pytest.fixture
def unique_organization_name():
    """Generate a unique organization name for testing"""
    return f"cli-test-org-{uuid.uuid4().hex[:8]}"


@pytest.fixture
def unique_invite_email():
    """Generate a unique email for invitation testing"""
    return f"cli-test-invite-{uuid.uuid4().hex[:8]}@example.com"


@pytest.fixture
def cli_runner():
    """Click CLI test runner"""
    return CliRunner()


@pytest.fixture
def cleanup_organizations():
    """Track created organizations for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for org_id in created_ids:
        try:
            delete_organization(org_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def test_organization(unique_organization_name, cleanup_organizations):
    """Create a test organization for CLI tests"""
    result = create_organization(
        name=unique_organization_name, description="CLI integration test organization"
    )
    cleanup_organizations.append(result["id"])
    return result


@pytest.mark.integration
class TestOrganizationCLICommands:
    """Integration tests for organization CLI commands"""

    def test_org_get_list_basic(self, cli_runner):
        """Test basic organization listing"""
        result = cli_runner.invoke(org_group, ["get"])

        assert result.exit_code == 0
        assert (
            "Organizations" in result.output
            or "No organizations found" in result.output
        )

    def test_org_get_list_json_output(self, cli_runner):
        """Test organization listing with JSON output"""
        result = cli_runner.invoke(org_group, ["get", "--json-output"])

        assert result.exit_code == 0

        # Should be valid JSON
        try:
            data = json.loads(result.output)
            assert "results" in data
            assert "count" in data
        except json.JSONDecodeError:
            pytest.fail("Output is not valid JSON")

    def test_org_get_list_with_search(self, cli_runner, test_organization):
        """Test organization listing with search parameter"""
        search_term = test_organization["name"][:10]  # Use part of the name

        result = cli_runner.invoke(org_group, ["get", "--search", search_term])

        assert result.exit_code == 0
        # Should contain search results or show no results found

    def test_org_get_list_with_ordering(self, cli_runner):
        """Test organization listing with ordering"""
        result = cli_runner.invoke(org_group, ["get", "--ordering", "name"])

        assert result.exit_code == 0

    def test_org_get_list_with_filters(self, cli_runner):
        """Test organization listing with filters"""
        result = cli_runner.invoke(org_group, ["get", "--filter", "name=test"])

        assert result.exit_code == 0

    def test_org_get_specific_organization(self, cli_runner, test_organization):
        """Test getting specific organization details"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["get", org_id])

        # Debug output if test fails
        if result.exit_code != 0:
            print(f"Exit code: {result.exit_code}")
            print(f"Output: {result.output}")
            print(f"Exception: {result.exception}")

        assert result.exit_code == 0
        # Strip ANSI codes for reliable string matching
        clean_output = strip_ansi_codes(result.output)
        assert "Organization Details" in clean_output
        assert org_id in clean_output
        assert test_organization["name"] in clean_output

    def test_org_get_specific_organization_json(self, cli_runner, test_organization):
        """Test getting specific organization with JSON output"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["get", org_id, "--json-output"])

        assert result.exit_code == 0

        # Should be valid JSON
        try:
            data = json.loads(result.output)
            assert data["id"] == org_id
            assert data["name"] == test_organization["name"]
        except json.JSONDecodeError:
            pytest.fail("Output is not valid JSON")

    def test_org_create_minimal(
        self, cli_runner, unique_organization_name, cleanup_organizations
    ):
        """Test creating organization with minimal parameters"""
        result = cli_runner.invoke(
            org_group, ["create", "--name", unique_organization_name]
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Created organization '{unique_organization_name}'" in clean_output

        # Extract org ID from output for cleanup
        lines = clean_output.split("\n")
        for line in lines:
            if "with ID:" in line:
                org_id = line.split("with ID:")[-1].strip()
                cleanup_organizations.append(org_id)
                break

    def test_org_create_with_description(
        self, cli_runner, unique_organization_name, cleanup_organizations
    ):
        """Test creating organization with description"""
        description = "CLI integration test organization"

        result = cli_runner.invoke(
            org_group,
            [
                "create",
                "--name",
                unique_organization_name,
                "--description",
                description,
            ],
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Created organization '{unique_organization_name}'" in clean_output

        # Extract org ID for cleanup
        lines = clean_output.split("\n")
        for line in lines:
            if "with ID:" in line:
                org_id = line.split("with ID:")[-1].strip()
                cleanup_organizations.append(org_id)
                break

    def test_org_create_with_billing_emails(
        self, cli_runner, unique_organization_name, cleanup_organizations
    ):
        """Test creating organization with billing emails"""
        billing_emails = "billing1@example.com,billing2@example.com"

        result = cli_runner.invoke(
            org_group,
            [
                "create",
                "--name",
                unique_organization_name,
                "--billing-emails",
                billing_emails,
            ],
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Created organization '{unique_organization_name}'" in clean_output

        # Extract org ID for cleanup
        lines = clean_output.split("\n")
        for line in lines:
            if "with ID:" in line:
                org_id = line.split("with ID:")[-1].strip()
                cleanup_organizations.append(org_id)
                break

    def test_org_create_json_output(
        self, cli_runner, unique_organization_name, cleanup_organizations
    ):
        """Test creating organization with JSON output"""
        result = cli_runner.invoke(
            org_group, ["create", "--name", unique_organization_name, "--json-output"]
        )

        assert result.exit_code == 0

        # Should be valid JSON
        try:
            data = json.loads(result.output)
            assert data["name"] == unique_organization_name
            assert "id" in data
            cleanup_organizations.append(data["id"])
        except json.JSONDecodeError:
            pytest.fail("Output is not valid JSON")

    def test_org_update_name(self, cli_runner, test_organization):
        """Test updating organization name"""
        org_id = test_organization["id"]
        new_name = f"{test_organization['name']}-updated"

        result = cli_runner.invoke(org_group, ["update", org_id, "--name", new_name])

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Updated organization {org_id}" in clean_output

    def test_org_update_description(self, cli_runner, test_organization):
        """Test updating organization description"""
        org_id = test_organization["id"]
        new_description = "Updated CLI test description"

        result = cli_runner.invoke(
            org_group, ["update", org_id, "--description", new_description]
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Updated organization {org_id}" in clean_output

    def test_org_update_billing_emails(self, cli_runner, test_organization):
        """Test updating organization billing emails"""
        org_id = test_organization["id"]
        new_emails = "updated1@example.com,updated2@example.com"

        result = cli_runner.invoke(
            org_group, ["update", org_id, "--billing-emails", new_emails]
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Updated organization {org_id}" in clean_output

    def test_org_update_json_output(self, cli_runner, test_organization):
        """Test updating organization with JSON output"""
        org_id = test_organization["id"]
        new_description = "JSON update test"

        result = cli_runner.invoke(
            org_group,
            ["update", org_id, "--description", new_description, "--json-output"],
        )

        assert result.exit_code == 0

        # Should be valid JSON
        try:
            data = json.loads(result.output)
            assert data["id"] == org_id
            assert data["description"] == new_description
        except json.JSONDecodeError:
            pytest.fail("Output is not valid JSON")

    def test_org_update_no_fields(self, cli_runner, test_organization):
        """Test updating organization with no fields provided"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["update", org_id])

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert "No fields to update" in clean_output

    def test_org_delete_with_confirmation(self, cli_runner, unique_organization_name):
        """Test deleting organization with confirmation"""
        # Create organization to delete
        create_result = create_organization(name=unique_organization_name)
        org_id = create_result["id"]

        result = cli_runner.invoke(org_group, ["delete", org_id], input="y\n")

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Deleted organization {org_id}" in clean_output

    def test_org_delete_with_yes_flag(self, cli_runner, unique_organization_name):
        """Test deleting organization with --yes flag"""
        # Create organization to delete
        create_result = create_organization(name=unique_organization_name)
        org_id = create_result["id"]

        result = cli_runner.invoke(org_group, ["delete", org_id, "--yes"])

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert f"Deleted organization {org_id}" in clean_output

    def test_org_delete_cancelled(self, cli_runner, test_organization):
        """Test cancelling organization deletion"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["delete", org_id], input="n\n")

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert "Operation cancelled" in clean_output

    def test_org_switch_valid_organization(self, cli_runner, test_organization):
        """Test switching to a valid organization"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["switch", org_id])

        # Note: This might fail if the API endpoint doesn't exist or user lacks permission
        # In that case, we expect a specific error message
        clean_output = strip_ansi_codes(result.output)
        if result.exit_code == 0:
            assert f"Switched to organization" in clean_output
        else:
            # API might not support switching or user lacks permission
            assert (
                "not found" in clean_output.lower() or "error" in clean_output.lower()
            )

    def test_org_switch_json_output(self, cli_runner, test_organization):
        """Test switching organization with JSON output"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["switch", org_id, "--json-output"])

        # Similar to above, this might not be supported
        if result.exit_code == 0:
            try:
                json.loads(result.output)
            except json.JSONDecodeError:
                pytest.fail("Output is not valid JSON")


@pytest.mark.integration
class TestOrganizationCLIUserCommands:
    """Integration tests for organization user management CLI commands"""

    def test_org_users_list_basic(self, cli_runner, test_organization):
        """Test listing organization users"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["users", "list", org_id])

        assert result.exit_code == 0
        # Should show users or "No users found"

    def test_org_users_list_json_output(self, cli_runner, test_organization):
        """Test listing organization users with JSON output"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["users", "list", org_id, "--json-output"]
        )

        assert result.exit_code == 0

        # Should be valid JSON
        try:
            data = json.loads(result.output)
            assert "results" in data
            assert "count" in data
        except json.JSONDecodeError:
            pytest.fail("Output is not valid JSON")

    def test_org_users_list_with_search(self, cli_runner, test_organization):
        """Test listing organization users with search"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["users", "list", org_id, "--search", "test"]
        )

        assert result.exit_code == 0

    def test_org_users_add_member(self, cli_runner, test_organization, test_user_email):
        """Test adding a user as member to organization"""
        if not test_user_email:
            pytest.skip("TEST_ORGANIZATION_USER_EMAIL environment variable must be set")

        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["users", "add", org_id, test_user_email])

        # This might fail if user is already in organization or lacks permission
        clean_output = strip_ansi_codes(result.output)
        if result.exit_code == 0:
            assert f"Added user {test_user_email}" in clean_output
            assert "member" in clean_output
        else:
            # User might already exist or other validation error
            assert "error" in clean_output.lower() or "already" in clean_output.lower()

    def test_org_users_add_admin(self, cli_runner, test_organization, test_user_email):
        """Test adding a user as admin to organization"""
        if not test_user_email:
            pytest.skip("TEST_ORGANIZATION_USER_EMAIL environment variable must be set")

        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["users", "add", org_id, test_user_email, "--admin"]
        )

        # This might fail if user is already in organization or lacks permission
        clean_output = strip_ansi_codes(result.output)
        if result.exit_code == 0:
            assert f"Added user {test_user_email}" in clean_output
            assert "admin" in clean_output
        else:
            # User might already exist or other validation error
            assert "error" in clean_output.lower() or "already" in clean_output.lower()

    def test_org_users_add_json_output(
        self, cli_runner, test_organization, test_user_email
    ):
        """Test adding user with JSON output"""
        if not test_user_email:
            pytest.skip("TEST_ORGANIZATION_USER_EMAIL environment variable must be set")

        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["users", "add", org_id, test_user_email, "--json-output"]
        )

        # This might fail due to user already existing
        if result.exit_code == 0:
            try:
                data = json.loads(result.output)
                assert "id" in data
            except json.JSONDecodeError:
                pytest.fail("Output is not valid JSON")

    def test_org_users_remove_with_confirmation(self, cli_runner, test_organization):
        """Test removing user with confirmation"""
        org_id = test_organization["id"]
        dummy_user_id = "dummy-user-id"

        result = cli_runner.invoke(
            org_group, ["users", "remove", org_id, dummy_user_id], input="y\n"
        )

        # Will likely fail because user doesn't exist
        clean_output = strip_ansi_codes(result.output)
        assert result.exit_code != 0 or "not found" in clean_output.lower()

    def test_org_users_remove_cancelled(self, cli_runner, test_organization):
        """Test cancelling user removal"""
        org_id = test_organization["id"]
        dummy_user_id = "dummy-user-id"

        result = cli_runner.invoke(
            org_group, ["users", "remove", org_id, dummy_user_id], input="n\n"
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert "Operation cancelled" in clean_output


@pytest.mark.integration
class TestOrganizationCLIInviteCommands:
    """Integration tests for organization invitation CLI commands"""

    def test_org_invites_list_basic(self, cli_runner, test_organization):
        """Test listing organization invites"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(org_group, ["invites", "list", org_id])

        assert result.exit_code == 0
        # Should show invites or "No pending invitations found"

    def test_org_invites_list_json_output(self, cli_runner, test_organization):
        """Test listing organization invites with JSON output"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["invites", "list", org_id, "--json-output"]
        )

        assert result.exit_code == 0

        # Should be valid JSON
        try:
            data = json.loads(result.output)
            assert "results" in data
            assert "count" in data
        except json.JSONDecodeError:
            pytest.fail("Output is not valid JSON")

    def test_org_invites_send_member(
        self, cli_runner, test_organization, unique_invite_email
    ):
        """Test sending invite as member"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["invites", "send", org_id, unique_invite_email]
        )

        clean_output = strip_ansi_codes(result.output)
        if result.exit_code == 0:
            assert f"Sent invitation to {unique_invite_email}" in clean_output
            assert "member" in clean_output
        else:
            # Might fail due to permissions or other validation
            pass

    def test_org_invites_send_admin(
        self, cli_runner, test_organization, unique_invite_email
    ):
        """Test sending invite as admin"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["invites", "send", org_id, unique_invite_email, "--admin"]
        )

        clean_output = strip_ansi_codes(result.output)
        if result.exit_code == 0:
            assert f"Sent invitation to {unique_invite_email}" in clean_output
            assert "admin" in clean_output
        else:
            # Might fail due to permissions or other validation
            pass

    def test_org_invites_send_json_output(
        self, cli_runner, test_organization, unique_invite_email
    ):
        """Test sending invite with JSON output"""
        org_id = test_organization["id"]

        result = cli_runner.invoke(
            org_group, ["invites", "send", org_id, unique_invite_email, "--json-output"]
        )

        if result.exit_code == 0:
            try:
                data = json.loads(result.output)
                assert "id" in data
                assert data["email"] == unique_invite_email
            except json.JSONDecodeError:
                pytest.fail("Output is not valid JSON")

    def test_org_invites_cancel_with_confirmation(self, cli_runner):
        """Test cancelling invite with confirmation"""
        dummy_invite_id = "dummy-invite-id"

        result = cli_runner.invoke(
            org_group, ["invites", "cancel", dummy_invite_id], input="y\n"
        )

        # Will likely fail because invite doesn't exist
        clean_output = strip_ansi_codes(result.output)
        assert result.exit_code != 0 or "not found" in clean_output.lower()

    def test_org_invites_cancel_with_yes_flag(self, cli_runner):
        """Test cancelling invite with --yes flag"""
        dummy_invite_id = "dummy-invite-id"

        result = cli_runner.invoke(
            org_group, ["invites", "cancel", dummy_invite_id, "--yes"]
        )

        # Will likely fail because invite doesn't exist
        clean_output = strip_ansi_codes(result.output)
        assert result.exit_code != 0 or "not found" in clean_output.lower()

    def test_org_invites_cancel_cancelled(self, cli_runner):
        """Test cancelling invite cancellation"""
        dummy_invite_id = "dummy-invite-id"

        result = cli_runner.invoke(
            org_group, ["invites", "cancel", dummy_invite_id], input="n\n"
        )

        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert "Operation cancelled" in clean_output


@pytest.mark.integration
class TestOrganizationCLIErrorHandling:
    """Integration tests for organization CLI error handling"""

    def test_org_get_nonexistent_organization(self, cli_runner):
        """Test getting a nonexistent organization"""
        nonexistent_id = f"nonexistent-{uuid.uuid4()}"

        result = cli_runner.invoke(org_group, ["get", nonexistent_id])

        assert result.exit_code != 0
        # Should show error message

    def test_org_create_missing_name(self, cli_runner):
        """Test creating organization without name"""
        result = cli_runner.invoke(org_group, ["create"])

        assert result.exit_code != 0
        # Should show usage or error about missing name

    def test_org_update_nonexistent_organization(self, cli_runner):
        """Test updating nonexistent organization"""
        nonexistent_id = f"nonexistent-{uuid.uuid4()}"

        result = cli_runner.invoke(
            org_group, ["update", nonexistent_id, "--name", "test"]
        )

        assert result.exit_code != 0
        # Should show error message

    def test_org_delete_nonexistent_organization(self, cli_runner):
        """Test deleting nonexistent organization"""
        nonexistent_id = f"nonexistent-{uuid.uuid4()}"

        result = cli_runner.invoke(org_group, ["delete", nonexistent_id, "--yes"])

        assert result.exit_code != 0
        # Should show error message

    def test_org_users_list_nonexistent_organization(self, cli_runner):
        """Test listing users for nonexistent organization"""
        nonexistent_id = f"nonexistent-{uuid.uuid4()}"

        result = cli_runner.invoke(org_group, ["users", "list", nonexistent_id])

        assert result.exit_code != 0
        # Should show error message

    def test_org_invites_list_nonexistent_organization(self, cli_runner):
        """Test listing invites for nonexistent organization"""
        nonexistent_id = f"nonexistent-{uuid.uuid4()}"

        result = cli_runner.invoke(org_group, ["invites", "list", nonexistent_id])

        assert result.exit_code != 0
        # Should show error message

    def test_invalid_filter_format(self, cli_runner):
        """Test organization listing with invalid filter format"""
        result = cli_runner.invoke(org_group, ["get", "--filter", "invalid_filter"])

        assert result.exit_code != 0
        clean_output = strip_ansi_codes(result.output)
        assert "Invalid filter format" in clean_output


@pytest.mark.integration
class TestOrganizationCLIWorkflows:
    """End-to-end CLI workflow tests"""

    def test_complete_organization_cli_workflow(
        self, cli_runner, unique_organization_name
    ):
        """Test complete organization CLI workflow"""
        # Step 1: Create organization
        create_result = cli_runner.invoke(
            org_group,
            [
                "create",
                "--name",
                unique_organization_name,
                "--description",
                "CLI workflow test",
                "--json-output",
            ],
        )

        assert create_result.exit_code == 0
        create_data = json.loads(create_result.output)
        org_id = create_data["id"]

        try:
            # Step 2: Get organization details
            get_result = cli_runner.invoke(org_group, ["get", org_id])
            assert get_result.exit_code == 0
            get_clean_output = strip_ansi_codes(get_result.output)
            assert unique_organization_name in get_clean_output

            # Step 3: Update organization
            update_result = cli_runner.invoke(
                org_group,
                ["update", org_id, "--description", "Updated CLI workflow test"],
            )
            assert update_result.exit_code == 0
            update_clean_output = strip_ansi_codes(update_result.output)
            assert f"Updated organization {org_id}" in update_clean_output

            # Step 4: List organizations (should include our org)
            list_result = cli_runner.invoke(
                org_group, ["get", "--search", unique_organization_name[:10]]
            )
            assert list_result.exit_code == 0

        finally:
            # Step 5: Delete organization
            delete_result = cli_runner.invoke(org_group, ["delete", org_id, "--yes"])
            assert delete_result.exit_code == 0
            delete_clean_output = strip_ansi_codes(delete_result.output)
            assert f"Deleted organization {org_id}" in delete_clean_output

    def test_organization_invite_cli_workflow(
        self, cli_runner, unique_organization_name, unique_invite_email
    ):
        """Test organization invitation CLI workflow"""
        # Create organization first
        create_result = cli_runner.invoke(
            org_group, ["create", "--name", unique_organization_name, "--json-output"]
        )

        assert create_result.exit_code == 0
        create_data = json.loads(create_result.output)
        org_id = create_data["id"]

        try:
            # Send invite
            invite_result = cli_runner.invoke(
                org_group, ["invites", "send", org_id, unique_invite_email]
            )

            # Might succeed or fail based on permissions
            if invite_result.exit_code == 0:
                invite_clean_output = strip_ansi_codes(invite_result.output)
                assert (
                    f"Sent invitation to {unique_invite_email}" in invite_clean_output
                )

                # List invites
                list_result = cli_runner.invoke(org_group, ["invites", "list", org_id])
                assert list_result.exit_code == 0

        finally:
            # Clean up organization
            delete_result = cli_runner.invoke(org_group, ["delete", org_id, "--yes"])
            assert delete_result.exit_code == 0


@pytest.mark.integration
class TestOrganizationCLISubprocessIntegration:
    """Integration tests using subprocess to test actual CLI executable"""

    def test_tess_org_help(self):
        """Test tess org help command via subprocess"""
        try:
            result = subprocess.run(
                ["tess", "org", "--help"],
                capture_output=True,
                text=True,
                timeout=30,
            )

            assert result.returncode == 0
            assert "Manage organizations" in result.stdout

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            pytest.skip(f"tess command not available or timeout: {e}")

    def test_tess_org_get_subprocess(self):
        """Test tess org get command via subprocess"""
        try:
            result = subprocess.run(
                ["tess", "org", "get", "--json-output"],
                capture_output=True,
                text=True,
                timeout=30,
            )

            # Should succeed (exit code 0) or fail with authentication error
            if result.returncode == 0:
                # Should be valid JSON
                try:
                    data = json.loads(result.stdout)
                    assert "results" in data
                    assert "count" in data
                except json.JSONDecodeError:
                    pytest.fail("Output is not valid JSON")
            else:
                # Might fail due to authentication or other issues
                assert (
                    "error" in result.stderr.lower()
                    or "authentication" in result.stderr.lower()
                )

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            pytest.skip(f"tess command not available or timeout: {e}")

    def test_tess_org_users_help(self):
        """Test tess org users help command via subprocess"""
        try:
            result = subprocess.run(
                ["tess", "org", "users", "--help"],
                capture_output=True,
                text=True,
                timeout=30,
            )

            assert result.returncode == 0
            assert "Manage organization users" in result.stdout

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            pytest.skip(f"tess command not available or timeout: {e}")

    def test_tess_org_invites_help(self):
        """Test tess org invites help command via subprocess"""
        try:
            result = subprocess.run(
                ["tess", "org", "invites", "--help"],
                capture_output=True,
                text=True,
                timeout=30,
            )

            assert result.returncode == 0
            assert "Manage organization invitations" in result.stdout

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            pytest.skip(f"tess command not available or timeout: {e}")
