import pytest
from atomict.io.cif import CIFAnalyzer
import os

# Test file path
TEST_FILE = os.path.join(os.path.dirname(__file__), './assets/GaAs.cif')

def test_cif_analyzer_initialization():
    analyzer = CIFAnalyzer(TEST_FILE)
    assert analyzer.atoms is not None
    assert analyzer.file_path == TEST_FILE

def test_basic_info():
    analyzer = CIFAnalyzer(TEST_FILE)
    basic_info = analyzer._get_basic_info()
    
    assert 'chemical_formula' in basic_info
    assert basic_info['chemical_formula'] in ['Ga4As4', 'As4Ga4']
    assert 'volume' in basic_info
    assert isinstance(basic_info['volume'], float)

def test_cell_parameters():
    analyzer = CIFAnalyzer(TEST_FILE)
    cell_info = analyzer._get_cell_parameters()
    
    assert 'cell_parameters' in cell_info
    params = cell_info['cell_parameters']
    
    # Check all parameters exist and are float
    for param in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
        assert param in params
        assert isinstance(params[param], float)
    
    # GaAs is cubic, so angles should be 90°
    assert pytest.approx(params['alpha'], 0.1) == 90.0
    assert pytest.approx(params['beta'], 0.1) == 90.0
    assert pytest.approx(params['gamma'], 0.1) == 90.0

def test_atomic_info():
    analyzer = CIFAnalyzer(TEST_FILE)
    atomic_info = analyzer._get_atomic_info()
    
    assert 'num_atoms' in atomic_info
    assert atomic_info['num_atoms'] == 8  # 4 Ga + 4 As
    
    assert 'element_counts' in atomic_info
    assert atomic_info['element_counts'] == {'Ga': 4, 'As': 4}
    
    assert 'atomic_numbers' in atomic_info
    assert len(atomic_info['atomic_numbers']) == 8
    
    assert 'positions' in atomic_info
    assert len(atomic_info['positions']) == 8

def test_symmetry_info():
    analyzer = CIFAnalyzer(TEST_FILE)
    symmetry_info = analyzer._get_symmetry_info()
    
    assert 'symmetry' in symmetry_info
    sym = symmetry_info['symmetry']
    
    # GaAs has F-43m space group (216)
    assert sym['space_group_number'] == 216
    assert sym['space_group_international'] == 'F-43m'
    
    # Should have 8 atoms in conventional cell
    assert sym['conventional_atoms'] == 8

def test_full_analysis():
    analyzer = CIFAnalyzer(TEST_FILE)
    analysis = analyzer.analyze()
    
    # Check that all main sections are present
    assert all(key in analysis for key in [
        'chemical_formula',
        'volume',
        'cell_parameters',
        'num_atoms',
        'element_counts',
        'symmetry'
    ])

def test_invalid_file():
    with pytest.raises(ValueError):
        CIFAnalyzer('nonexistent.cif')
