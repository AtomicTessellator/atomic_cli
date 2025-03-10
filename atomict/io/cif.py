import logging
import ase.io
import spglib
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure


class CIFAnalyzer:
    def __init__(self, file_path):
        self.file_path = file_path
        self.atoms = None
        self._load_structure()

    def _load_structure(self):
        try:
            self.atoms = ase.io.read(self.file_path, format='cif')
            if self.atoms is None:
                raise ValueError('Could not parse CIF file structure')
        except Exception as e:
            raise ValueError(f'Not a valid CIF file: {str(e)}')

    def analyze(self):
        """Performs complete analysis of the CIF structure"""
        analysis = {}
        
        # Basic structure info
        analysis.update(self._get_basic_info())
        
        # Cell parameters
        analysis.update(self._get_cell_parameters())
        
        # Atomic information
        analysis.update(self._get_atomic_info())
        
        # Symmetry analysis
        analysis.update(self._get_symmetry_info())
        
        return analysis

    def _get_basic_info(self):
        info = {}
        try:
            info['chemical_formula'] = self.atoms.get_chemical_formula()
        except Exception as e:
            info['chemical_formula'] = 'Unknown'
            logging.warning(f'Could not determine chemical formula: {str(e)}')
        
        try:
            info['volume'] = float(self.atoms.get_volume())
        except Exception:
            info['volume'] = None
            
        return info

    def _get_cell_parameters(self):
        try:
            cell_angles = self.atoms.cell.angles()
            return {
                'cell_parameters': {
                    'a': float(self.atoms.cell[0][0]) if self.atoms.cell[0][0] else None,
                    'b': float(self.atoms.cell[1][1]) if self.atoms.cell[1][1] else None,
                    'c': float(self.atoms.cell[2][2]) if self.atoms.cell[2][2] else None,
                    'alpha': float(cell_angles[0]) if cell_angles[0] else None,
                    'beta': float(cell_angles[1]) if cell_angles[1] else None,
                    'gamma': float(cell_angles[2]) if cell_angles[2] else None,
                }
            }
        except Exception as e:
            logging.warning(f'Could not determine cell parameters: {str(e)}')
            return {'cell_parameters': None}

    def _get_atomic_info(self):
        try:
            info = {
                'num_atoms': len(self.atoms),
                'atomic_numbers': self.atoms.numbers.tolist(),
                'chemical_symbols': self.atoms.get_chemical_symbols(),
                'positions': self.atoms.positions.tolist(),
                'cell': self.atoms.cell.tolist(),
                'pbc': self.atoms.pbc.tolist()
            }
            
            # Get unique elements and their counts
            unique_elements = {}
            for symbol in self.atoms.get_chemical_symbols():
                unique_elements[symbol] = unique_elements.get(symbol, 0) + 1
            info['element_counts'] = unique_elements
            
            return info
        except Exception as e:
            logging.warning(f'Could not get complete atomic information: {str(e)}')
            return {
                'num_atoms': 0,
                'atomic_numbers': [],
                'chemical_symbols': [],
                'positions': [],
                'cell': None,
                'pbc': None,
                'element_counts': {}
            }

    def _get_symmetry_info(self):
        try:
            cell = (self.atoms.cell, self.atoms.get_scaled_positions(), self.atoms.numbers)
            dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5)
            primitive_cell = spglib.find_primitive(cell, symprec=1e-5)
            
            if dataset:
                unique_equiv = set(dataset.equivalent_atoms)
                wyckoff_multiplicities = {
                    letter: list(dataset.wyckoffs).count(letter)
                    for letter in set(dataset.wyckoffs)
                }
                
                # Convert ASE atoms to pymatgen Structure for SpacegroupAnalyzer
                lattice = self.atoms.get_cell()
                species = self.atoms.get_chemical_symbols()
                coords = self.atoms.get_scaled_positions()
                structure = Structure(lattice, species, coords)
                
                # Now use SpacegroupAnalyzer with the pymatgen Structure
                analyzer = SpacegroupAnalyzer(structure)
                bravais = analyzer.get_lattice_type()

                symmetry_info = {
                    'symmetry': {
                        'space_group_number': dataset.number,
                        'space_group_international': dataset.international,
                        'hall_symbol': dataset.hall,
                        'primitive_atoms': len(primitive_cell[2]) if primitive_cell else len(self.atoms),
                        'conventional_atoms': len(self.atoms),
                        'wyckoff_letters': dataset.wyckoffs,
                        'equivalent_atoms': dataset.equivalent_atoms.tolist(),
                        'transformation_matrix': dataset.transformation_matrix.tolist(),
                        'origin_shift': dataset.origin_shift.tolist(),
                        'asymmetric_unit_atoms': len(unique_equiv),
                        'wyckoff_multiplicities': wyckoff_multiplicities,
                        'bravais_lattice': bravais
                    }
                }
                return symmetry_info
        except Exception as e:
            logging.warning(f'Could not determine symmetry information: {str(e)}')
        
        return {'symmetry': None}
