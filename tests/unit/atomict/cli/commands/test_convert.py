import os
import pytest
from unittest.mock import patch, MagicMock
from click.testing import CliRunner
from pathlib import Path
from atomict.cli.main import convert_command


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def fixture_path():
    """Return the path to the fixtures directory."""
    return Path(__file__).parent.parent.parent.parent.parent / "fixtures"


@pytest.fixture
def mock_ase_modules():
    # Create a comprehensive mock for ASE and related modules
    with patch.dict('sys.modules', {
        'ase': MagicMock(),
        'ase.io': MagicMock(),
        'ase.io.read': MagicMock(),
        'ase.io.write': MagicMock(),
        'ase.io.formats': MagicMock(),
        'atomict.io.msgpack': MagicMock()
    }) as mocked_modules:
        # Configure the atom object that will be returned by read
        mock_atoms = MagicMock()
        mocked_modules['ase.io'].read.return_value = mock_atoms
        
        # Configure the UnknownFileTypeError exception
        mocked_modules['ase.io.formats'].UnknownFileTypeError = type(
            'UnknownFileTypeError', (Exception,), {})
        
        yield mocked_modules


@pytest.mark.parametrize("input_fixture,output_format", [
    ("Al.cif", "xyz"),
    ("Be.cif", "xyz"),
    ("Al.cif", "atm"),
])
def test_successful_conversion(runner, mock_ase_modules, input_fixture, output_format, fixture_path, tmp_path):
    # Use real fixture files
    input_file = fixture_path / input_fixture
    output_file = tmp_path / f"output.{output_format}"
    
    # Make sure the fixture exists
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
        
    with patch('os.path.exists', return_value=True):
        result = runner.invoke(convert_command, [str(input_file), str(output_file)])
        
    # Check that the command ran successfully
    assert result.exit_code == 0
    assert "Successfully converted" in result.output
    
    # Verify the appropriate functions were called
    mock_ase_modules['ase.io'].read.assert_called_once_with(str(input_file))
    
    # If output is .atm format, should use msgpack
    if output_format == 'atm':
        mock_ase_modules['atomict.io.msgpack'].save_msgpack.assert_called_once()
    else:
        mock_ase_modules['ase.io'].write.assert_called_once()


def test_input_file_not_found(runner, mock_ase_modules, tmp_path):
    nonexistent_file = tmp_path / "nonexistent.xyz"
    output_file = tmp_path / "output.cif"
    
    with patch('os.path.exists', return_value=False):
        result = runner.invoke(convert_command, [str(nonexistent_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "Input file" in result.output
    assert str(nonexistent_file) in result.output
    
    # Read should not be called if the file doesn't exist
    mock_ase_modules['ase.io'].read.assert_not_called()


def test_unsupported_input_format(runner, mock_ase_modules, tmp_path):
    # Create a test file with unsupported extension
    unsupported_file = tmp_path / "test.unsupported"
    unsupported_file.write_text("dummy content")
    output_file = tmp_path / "output.cif"
    
    with patch('os.path.exists', return_value=True):
        result = runner.invoke(convert_command, [str(unsupported_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "not supported" in result.output
    
    # Read should not be called if the format is unsupported
    mock_ase_modules['ase.io'].read.assert_not_called()


def test_unsupported_output_format(runner, mock_ase_modules, fixture_path, tmp_path):
    # Use a real fixture file
    input_file = fixture_path / "Al.cif"
    output_file = tmp_path / "output.unsupported"
    
    # Skip if fixture doesn't exist
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
        
    with patch('os.path.exists', return_value=True):
        result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "not supported" in result.output
    
    # Read might be called but write should not
    mock_ase_modules['ase.io'].write.assert_not_called()


def test_read_error(runner, mock_ase_modules, fixture_path, tmp_path):
    # Use a real fixture file
    input_file = fixture_path / "Al.cif"
    output_file = tmp_path / "output.xyz"
    
    # Skip if fixture doesn't exist
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
        
    # Configure read to raise an exception
    unknown_file_type_error = mock_ase_modules['ase.io.formats'].UnknownFileTypeError
    mock_ase_modules['ase.io'].read.side_effect = unknown_file_type_error("Test error")
    
    with patch('os.path.exists', return_value=True):
        result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "Unknown file type" in result.output


def test_write_error(runner, mock_ase_modules, fixture_path, tmp_path):
    # Use a real fixture file
    input_file = fixture_path / "Al.cif"
    output_file = tmp_path / "output.xyz"
    
    # Skip if fixture doesn't exist
    if not input_file.exists():
        pytest.skip(f"Fixture file {input_file} not found")
        
    # Configure write to raise an exception
    mock_ase_modules['ase.io'].write.side_effect = Exception("Write error")
    
    with patch('os.path.exists', return_value=True):
        result = runner.invoke(convert_command, [str(input_file), str(output_file)])
    
    assert result.exit_code == 0  # The command exits gracefully
    assert "Error writing" in result.output
