import pytest
from click.testing import CliRunner
from pathlib import Path
# writes to tempfile due to conflicts with library code
import tempfile
from atomict.cli.main import convert_command


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def fixture_path():
    """Return the path to the fixtures directory."""
    return Path(__file__).resolve().parents[4] / "fixtures"


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)


@pytest.mark.parametrize("input_fixture,output_format", [
    ("Al.cif", "xyz"),
    ("Be.cif", "xyz"),
    ("Al.cif", "atm"),
])
def test_successful_conversion(runner, fixture_path, temp_dir, input_fixture, output_format):
    # Use real fixture files
    input_file = fixture_path / input_fixture
    output_file = temp_dir / f"output.{output_format}"
    
    # Make sure the fixture exists
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
    
    # Run the conversion command
    result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    # Check that the command ran successfully
    assert result.exit_code == 0
    assert "Successfully converted" in result.output
    
    # Verify the output file was created
    assert output_file.exists()


def test_input_file_not_found(runner, temp_dir):
    nonexistent_file = temp_dir / "nonexistent.xyz"
    output_file = temp_dir / "output.cif"
    
    result = runner.invoke(convert_command, [str(nonexistent_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "Input file" in result.output
    assert str(nonexistent_file) in result.output
    
    # Verify output file was not created
    assert not output_file.exists()


def test_unsupported_input_format(runner, temp_dir):
    # Create a test file with unsupported extension
    unsupported_file = temp_dir / "test.unsupported"
    unsupported_file.write_text("dummy content")
    output_file = temp_dir / "output.cif"
    
    result = runner.invoke(convert_command, [str(unsupported_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "not supported" in result.output
    
    # Verify output file was not created
    assert not output_file.exists()


def test_unsupported_output_format(runner, fixture_path, temp_dir):
    # Use a real fixture file
    input_file = fixture_path / "Al.cif"
    output_file = temp_dir / "output.unsupported"
    
    # Skip if fixture doesn't exist
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
    
    result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "not supported" in result.output
    
    # Verify output file was not created
    assert not output_file.exists()


def test_read_error(runner, fixture_path, temp_dir, monkeypatch):
    # Use a real fixture file
    input_file = fixture_path / "Al.cif"
    output_file = temp_dir / "output.xyz"
    
    # Skip if fixture doesn't exist
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
    
    # Create a mock function that raises an exception
    def mock_read(*args, **kwargs):
        raise Exception("Unknown file type")
    
    # Patch the read function to raise an exception
    import ase.io
    monkeypatch.setattr(ase.io, "read", mock_read)
    
    result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "Unknown file type" in result.output
    
    # Verify output file was not created
    assert not output_file.exists()


def test_write_error(runner, fixture_path, temp_dir, monkeypatch):
    # Use a real fixture file
    input_file = fixture_path / "Al.cif"
    output_file = temp_dir / "output.xyz"
    
    # Skip if fixture doesn't exist
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
    
    # Create a mock atoms object and read function that returns it
    import ase
    mock_atoms = object()  # Simple mock atoms object
    
    def mock_read(*args, **kwargs):
        return mock_atoms
    
    # Create a mock write function that raises an exception
    def mock_write(*args, **kwargs):
        raise Exception("Write error")
    
    # Patch the functions
    import ase.io
    monkeypatch.setattr(ase.io, "read", mock_read)
    monkeypatch.setattr(ase.io, "write", mock_write)
    
    result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "Error writing" in result.output
    
    # Verify output file was not created
    assert not output_file.exists()
