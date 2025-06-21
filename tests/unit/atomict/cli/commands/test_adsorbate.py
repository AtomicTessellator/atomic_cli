from unittest.mock import Mock, patch

import pytest
from click.testing import CliRunner

from atomict.cli.commands.adsorbate import create, delete, get, update


@pytest.fixture
def mock_client():
    with patch("atomict.cli.commands.adsorbate.get_client") as mock:
        client = Mock()
        mock.return_value = client
        yield client


@pytest.fixture
def runner():
    return CliRunner()


def test_get_single_adsorbate(runner, mock_client):
    # Mock response data
    mock_client.get.return_value = {
        "id": "123",
        "smiles": "CC(=O)O",
        "binding_indices": [1, 2],
    }

    result = runner.invoke(get, ["123"])
    assert result.exit_code == 0
    assert "123" in result.output
    mock_client.get.assert_called_once_with("/api/adsorbate/123/")


def test_get_all_adsorbates(runner, mock_client):
    # Mock response data
    mock_client.get_all.return_value = [
        {"id": "123", "smiles": "CC(=O)O"},
        {"id": "456", "smiles": "CH4"},
    ]

    result = runner.invoke(get, ["--all"])

    assert result.exit_code == 0
    mock_client.get_all.assert_called_once_with("/api/adsorbate/", params={})


def test_create_adsorbate(runner, mock_client):
    mock_client.post.return_value = {"id": "789"}

    # Either smiles OR ase_atoms is required, not both
    result = runner.invoke(
        create,
        [
            "--ase-atoms",
            "Atoms(...)",  # Changed to use ase_atoms
            "--binding-indices",
            "1",
            "--binding-indices",
            "2",
        ],
    )

    assert result.exit_code == 0
    assert "789" in result.output


def test_update_adsorbate(runner, mock_client):
    # Mock response data
    mock_client.put.return_value = {"id": "123"}

    result = runner.invoke(
        update, ["123", "--smiles", "CH4", "--reaction-string", "A -> B"]
    )

    assert result.exit_code == 0
    assert "123" in result.output
    mock_client.put.assert_called_once_with(
        "/api/adsorbate/123/", {"smiles": "CH4", "reaction_string": "A -> B"}
    )


def test_delete_adsorbate(runner, mock_client):
    result = runner.invoke(delete, ["123"])

    assert result.exit_code == 0
    assert "123" in result.output
    mock_client.delete.assert_called_once_with("/api/adsorbate/123/")


def test_get_json_output(runner, mock_client):
    mock_client.get.return_value = {"id": "123", "smiles": "CC(=O)O"}

    result = runner.invoke(get, ["123", "--json-output"])
    assert result.exit_code == 0
    assert '"id": "123"' in result.output
    assert '"smiles": "CC(=O)O"' in result.output


def test_create_adsorbate_with_ase_atoms(runner, mock_client):
    mock_client.post.return_value = {"id": "789"}

    result = runner.invoke(
        create, ["--ase-atoms", "Atoms(...)", "--reaction-string", "A + B -> C"]
    )

    assert result.exit_code == 0
    assert "789" in result.output
    mock_client.post.assert_called_once_with(
        "/api/adsorbate/", {"ase_atoms": "Atoms(...)", "reaction_string": "A + B -> C"}
    )


def test_create_validates_binding_indices(runner, mock_client):
    mock_client.post.return_value = {"id": "789"}

    result = runner.invoke(
        create,
        [
            "--ase-atoms",
            "Atoms(...)",  # Added required ase_atoms
            "--binding-indices",
            "1",
            "--binding-indices",
            "2",
        ],
    )

    assert result.exit_code == 0
    mock_client.post.assert_called_with(
        "/api/adsorbate/", {"ase_atoms": "Atoms(...)", "binding_indices": [1, 2]}
    )


def test_create_handles_missing_required_args(runner, mock_client):
    """Test that proper error is shown when required args are missing"""
    result = runner.invoke(create)

    assert result.exit_code != 0
    assert "Error: Missing option" in result.output  # More generic assertion


def test_get_formats_json_output(runner, mock_client):
    """Test that --json-output flag changes the output format"""
    mock_client.get.return_value = {"id": "123", "smiles": "CC(=O)O"}

    # Test normal output
    result = runner.invoke(get, ["123"])
    assert "ID: 123" in result.output

    # Test JSON output
    json_result = runner.invoke(get, ["123", "--json-output"])
    assert '"id": "123"' in json_result.output
