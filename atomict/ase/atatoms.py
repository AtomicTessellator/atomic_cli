import uuid
import time
import numpy as np
import json
import datetime
import requests
import os
import sys
import functools
import hashlib
from typing import Optional, Dict, Any, List, Union, Tuple
import logging

from ase.atoms import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT

from deepdiff import DeepDiff
from atomict.api import post


logger = logging.getLogger(__name__)
# Configure logger to print to stdout
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

def track_state_changes():
    """
    Decorator to track state changes after method execution.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            # Call the original method
            result = func(self, *args, **kwargs)
            logger.info(f"got {func.__name__} result: {result}")
            # Track state changes after the method call
            if hasattr(self, '_capture_state_diff'):
                diff = self._capture_state_diff()
                # Send diff right after capturing it
                if hasattr(self, '_send_diff') and diff:
                    logger.info(f"posting diff: {diff}")
                    self._send_diff(diff)
            
            return result
        return wrapper
    return decorator

class ATAtoms:
    """
    A lightweight wrapper around ASE Atoms that transparently tracks modifications
    and generates diffs for server synchronization.
    """
    
    def __init__(self, atoms: Atoms, server_url: Optional[str] = None, 
                 batch_size: int = 10, sync_interval: float = 60.0):
        """
        Initialize the ATAtoms wrapper.
        
        Parameters:
        -----------
        atoms: ASE Atoms object to wrap
        server_url: URL of the server to send diffs to (if None, will use AT_SERVER env var)
        batch_size: Number of diffs to accumulate before sending to server
        sync_interval: Maximum time in seconds between syncs to server
        """
        # Validate that atoms is an ASE Atoms object
        if not isinstance(atoms, Atoms):
            raise TypeError(f"Expected ASE Atoms object, got {type(atoms).__name__}: {atoms}")
            
        # Initialize core properties
        self._atoms = atoms
        self._server_url = server_url or os.environ.get('AT_SERVER')
        self._object_id = str(uuid.uuid4())
        self._batch_size = batch_size
        self._sync_interval = sync_interval
        self._diffs = []  # Track the state diffs
        self._last_sync_time = time.time()
        self._seq_num = 0  # Initialize sequence number counter
        self._state_id = None  # Will be set after server initialization
        self._initialized = False  # Flag to track initialization status
        
        # Initialize state tracking
        self._initial_state = self._get_current_state()  # Save the initial state
        self._previous_state = self._serialize_state(self._initial_state)  # Serialize for comparison
        
        # Generate structure_id hash
        self._structure_id = self._hash_state(self._previous_state)
        
        # Record the initial state
        self._capture_state_diff()
    
    def _hash_state(self, state_data):
        """Generate a SHA-256 hash of the state data"""
        # Sort the keys to ensure consistent hash for same data
        state_json = json.dumps(state_data, sort_keys=True)
        _hash = hashlib.sha256(state_json.encode('utf-8')).hexdigest()
        logger.info(f"created structure id: {_hash}")
        return _hash
    
    def _initialize_on_server(self):
        """Send the initial state to create an object on the server"""
        # Try to use self._server_url or fall back to environment variable
        if not self._server_url:
            self._server_url = os.environ.get('AT_SERVER')
            if not self._server_url:
                print("Warning: No server URL provided and AT_SERVER environment variable not set")
                return
        
        try:
            # Extract and format the initial state according to server model
            state_data = {
                'id': self._object_id,
                'structure_id': self._structure_id,
                'numbers': self._previous_state.get('numbers'),
                'positions': self._previous_state.get('positions'),
                'cell': self._previous_state.get('cell'),
                'pbc': self._previous_state.get('pbc'),
                'energy': self._previous_state.get('energy'),
                'symbols': self._previous_state.get('symbols'),
                'forces': self._previous_state.get('forces'),
                'stress': self._previous_state.get('stress'),
                'info': self._previous_state.get('info', {}),
                'scaled_positions': self._previous_state.get('scaled_positions'),
                'tags': self._previous_state.get('tags'),
                'momenta': self._previous_state.get('momenta'),
                'velocities': self._previous_state.get('velocities'),
                'masses': self._previous_state.get('masses'),
                'magmoms': self._previous_state.get('magmoms'),
                'charges': self._previous_state.get('charges'),
                'celldisp': self._previous_state.get('celldisp'),
                'constraints': self._previous_state.get('constraints')
            }
            # Use the API module to make the request
            response = post(
                'api/atatoms-states', 
                state_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            # Store the state ID for future reference
            if response and 'id' in response:
                self._state_id = response['id']
                print(f"Successfully initialized state on server with ID: {self._state_id}")
            else:
                print(f"Warning: Server response did not contain expected 'id' field: {response}")
            
        except Exception as e:
            print(f"Warning: Failed to initialize object on server: {e}")
    
    def _capture_state_diff(self):
        """
        Capture differences between current and previous state and append to diffs array
        """
        current_state = self._get_current_state()
        serialized_current = self._serialize_state(current_state)
        
        # For initial state, just initialize on server without creating a diff
        if not self._initialized:
            # Initialize on server first to get state_id
            self._initialize_on_server()
            # Mark as initialized
            self._initialized = True
            return None  # No diff for initial state
        
        # Generate diff using DeepDiff
        diff = DeepDiff(self._previous_state, serialized_current, verbose_level=1)
        
        # Only record if there are actual changes
        if diff:
            timestamp = datetime.datetime.now().isoformat()
            
            # Record the diff for state tracking
            self._diffs.append({
                'timestamp': timestamp,
                'sequence': self._seq_num,
                'diff': diff
            })
            
            # Update previous state for next comparison
            self._previous_state = serialized_current
            self._seq_num += 1  # Increment sequence number
            
            # Send diff only if we have a state_id from initialization
            if self._state_id:
                return diff
            else:
                print("Warning: Cannot send diff before initializing state on server")
                # Try to initialize again
                self._initialize_on_server()
                # If initialization succeeds, return diff to be sent
                if self._state_id:
                    return diff
        return None
    
    # Add special property for calculator
    @property
    def calc(self):
        return self._atoms.calc
    
    @calc.setter
    def calc(self, calculator):
        self._atoms.calc = calculator
        self._capture_state_diff()
    
    def _get_current_state(self) -> Dict[str, Any]:
        """Get the complete state of the atoms object"""
        atoms = self._atoms
        
        # Check if atoms is actually an ASE Atoms object, not a string or other type
        if isinstance(atoms, str):
            raise TypeError(f"Expected ASE Atoms object, got string: {atoms}")
        
        # Use atoms.todict() if available
        if hasattr(atoms, 'todict'):
            state = atoms.todict()
            # Ensure symbols are always included even if todict() doesn't include them
            if 'symbols' not in state:
                state['symbols'] = atoms.get_chemical_symbols()
        else:
            # Fallback to manual dictionary creation
            state = {}
            
            # Handle positions
            if hasattr(atoms, 'positions') and atoms.positions is not None:
                state['positions'] = atoms.positions.copy()
                
            # Handle cell
            if hasattr(atoms, 'cell'):
                if hasattr(atoms.cell, 'array'):
                    state['cell'] = atoms.cell.array.copy()
                else:
                    state['cell'] = atoms.cell.copy() if hasattr(atoms.cell, 'copy') else atoms.cell
                
            # Handle pbc
            if hasattr(atoms, 'pbc'):
                state['pbc'] = atoms.pbc.copy() if hasattr(atoms.pbc, 'copy') else atoms.pbc
            
            # Handle atomic numbers and symbols
            if hasattr(atoms, 'numbers'):
                state['numbers'] = atoms.numbers.copy()
                
            if hasattr(atoms, 'get_chemical_symbols'):
                state['symbols'] = atoms.get_chemical_symbols()
        
        # Add calculated properties if available
        if hasattr(atoms, 'get_potential_energy'):
            try:
                state['energy'] = atoms.get_potential_energy()
            except:
                pass
                
        if hasattr(atoms, 'get_forces'):
            try:
                state['forces'] = atoms.get_forces()
            except:
                pass
                
        if hasattr(atoms, 'get_stress'):
            try:
                state['stress'] = atoms.get_stress()
            except:
                pass
        
        # Add info dictionary to align with backend model
        if hasattr(atoms, 'info') and atoms.info:
            state['info'] = atoms.info.copy()
        
        # Add any custom arrays that might have been set
        for name in atoms.arrays:
            if name not in ['positions', 'numbers']:
                state[name] = atoms.arrays[name].copy()
                
        return state

    def _serialize_state(self, state: Dict[str, Any]) -> Dict[str, Any]:
        """Convert numpy arrays to lists for serialization and comparison"""
        serialized = {}
        for key, value in state.items():
            if isinstance(value, np.ndarray):
                serialized[key] = value.tolist()
            elif isinstance(value, (np.integer, np.floating)):
                serialized[key] = float(value) if isinstance(value, np.floating) else int(value)
            else:
                serialized[key] = value
        
        return serialized
    
    def _serialize_diff(self, diff):
        """Convert DeepDiff output to JSON serializable format"""
        if not diff:
            return diff
            
        # Create a new dict for the serialized diff
        serialized_diff = {}
        
        # Process each diff type
        for diff_type, values in diff.items():
            # Handle SetOrdered and other special types
            if hasattr(values, 'items'):  # For dictionary-like diff types
                serialized_diff[diff_type] = dict(values)
            elif isinstance(values, (list, tuple)):  # For lists
                serialized_diff[diff_type] = list(values)
            elif hasattr(values, '__iter__') and not isinstance(values, str):  # For other iterables
                serialized_diff[diff_type] = list(values)
            else:
                serialized_diff[diff_type] = values
                
        return serialized_diff
    
    # Special ASE methods that need explicit tracking
    @track_state_changes()
    def rotate(self, *args, **kwargs):
        result = self._atoms.rotate(*args, **kwargs)
        return self if result is self._atoms else result
    
    @track_state_changes()
    def translate(self, displacement):
        result = self._atoms.translate(displacement)
        return self if result is self._atoms else result
    
    @track_state_changes()
    def repeat(self, *args, **kwargs):
        result = self._atoms.repeat(*args, **kwargs)
        if result is self._atoms:
            return self
        else:
            wrapped = ATAtoms(result)
            return wrapped
            
    # Explicitly track rattle which modifies positions directly
    @track_state_changes()
    def rattle(self, *args, **kwargs):
        result = self._atoms.rattle(*args, **kwargs)
        return self if result is self._atoms else result
    
    # Handle individual atom access with change tracking
    def __getitem__(self, index):
        if isinstance(index, slice):
            # Handle slices by returning a new wrapped atoms
            new_atoms = self._atoms[index]
            wrapped = ATAtoms(new_atoms)
            return wrapped
        else:
            # Return a proxy for individual atom access
            return _AtomProxy(self, index)
    
    def __setitem__(self, index, value):
        self._atoms[index] = value
        self._capture_state_diff()
    
    def __delitem__(self, index):
        del self._atoms[index]
        self._capture_state_diff()
    
    def __iadd__(self, other):
        self._atoms += other
        self._capture_state_diff()
        return self
    
    def __len__(self):
        return len(self._atoms)
    
    # Transparent method forwarding
    def __getattr__(self, name):
        attr = getattr(self._atoms, name)
        
        if callable(attr):
            def wrapped_method(*args, **kwargs):
                result = attr(*args, **kwargs)
                # Handle methods that return the atoms object
                if result is self._atoms:
                    # Capture state changes after method execution
                    self._capture_state_diff()
                    return self
                elif isinstance(result, Atoms):
                    return ATAtoms(result)
                else:
                    return result
                    
            return wrapped_method
        else:
            return attr
    
    def __dir__(self):
        return list(set(dir(self.__class__)) | set(dir(self._atoms)))
    
    def __del__(self):
        """Ensure any remaining diffs are synced before garbage collection"""
        if hasattr(self, '_diffs') and self._diffs:
            try:
                self._sync_diffs()
            except:
                pass

    def _send_diff(self, diff):
        """
        Send a diff to the server at /api/atatoms-diffs
        """
        # Don't attempt to send if we don't have a state_id yet
        if not self._state_id:
            print("Warning: Cannot send diff without state_id. Attempting to initialize on server first.")
            self._initialize_on_server()
            if not self._state_id:
                return None
        
        try:
            # Get current state to include in the request
            current_state = self._get_current_state()
            serialized_current = self._serialize_state(current_state)
            
            # Serialize the diff to make it JSON compatible
            serialized_diff = self._serialize_diff(diff)
            
            # Prepare the data to send
            diff_data = {
                'atoms_id': self._object_id,
                'structure_id': self._structure_id,
                'state_id': self._state_id,
                'sequence': self._seq_num - 1,  # Use the sequence number assigned to this diff
                'timestamp': datetime.datetime.now().isoformat(),
                'data': serialized_diff,
                'state': self._state_id  # Use state_id as the reference to the state
            }
            
            # Send the diff to the server
            response = post(
                'api/atatoms-diffs',
                diff_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            if not response:
                print("Warning: Empty response when sending diff to server")
            
            return response
        except Exception as e:
            print(f"Warning: Failed to send diff to server: {e}")
            return None

    def save_current_state(self) -> Dict[str, Any]:
        """
        Get current state, serialize it, hash it, and send it to server.
        
        Returns:
        --------
        Dict with server response data
        """
        # Get and serialize current state
        current_state = self._get_current_state()
        serialized_state = self._serialize_state(current_state)
        
        # Hash the state
        structure_id = self._hash_state(serialized_state)
        
        # Prepare the data to send to server
        state_data = {
            'structure_id': structure_id,
            'numbers': serialized_state.get('numbers'),
            'positions': serialized_state.get('positions'),
            'cell': serialized_state.get('cell'),
            'pbc': serialized_state.get('pbc'),
            'energy': serialized_state.get('energy'),
            'symbols': serialized_state.get('symbols'),
            'forces': serialized_state.get('forces'),
            'stress': serialized_state.get('stress'),
            'info': serialized_state.get('info', {}),
            'scaled_positions': serialized_state.get('scaled_positions'),
            'tags': serialized_state.get('tags'),
            'momenta': serialized_state.get('momenta'),
            'velocities': serialized_state.get('velocities'),
            'masses': serialized_state.get('masses'),
            'magmoms': serialized_state.get('magmoms'),
            'charges': serialized_state.get('charges'),
            'celldisp': serialized_state.get('celldisp'),
            'constraints': serialized_state.get('constraints')
        }
        
        try:
            # Send the state to server
            logger.info(f"Saving current state with structure_id: {structure_id}")
            response = post(
                'api/atatoms-states', 
                state_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            if response and 'id' in response:
                logger.info(f"Successfully saved state on server with ID: {response['id']}")
                return response
            else:
                logger.warning(f"Server response did not contain expected 'id' field: {response}")
                return response
                
        except Exception as e:
            logger.error(f"Failed to save state on server: {e}")
            return {"error": str(e)}

    @classmethod
    def from_known_state(cls, state_id: str) -> 'ATAtoms':
        """Retrieves a state from the server and initializes an ATAtoms object with it."""
        raise NotImplementedError

    def track_changes(self):
        """
        Simple callback for ASE optimizers to track changes.
        Calculates diff between previous and current state and sends if changes exist.
        
        Usage:
            opt = BFGS(atoms)
            opt.attach(atoms.track_changes, interval=1)
            opt.run(fmax=0.02)
        """
        diff = self._capture_state_diff()
        
        if diff and hasattr(self, '_send_diff'):
            logger.info(f"posting diff from optimization step")
            self._send_diff(diff)
            
        return True


class _AtomProxy:
    """Proxy for individual atom access that tracks position changes"""
    
    def __init__(self, parent, index):
        self._parent = parent
        self._index = index
    
    @property
    def position(self):
        # Return a position proxy that tracks changes
        return _PositionProxy(self._parent, self._index)
    
    @position.setter
    def position(self, value):
        self._parent._atoms.positions[self._index] = value
        diff = self._parent._capture_state_diff()
        # Send diff if it exists
        if diff and hasattr(self._parent, '_send_diff'):
            logger.info(f"posting diff from atom proxy: {diff}")
            self._parent._send_diff(diff)
    
    def __getattr__(self, name):
        # Forward all other attributes to the actual atom
        return getattr(self._parent._atoms[self._index], name)
    
    def __setattr__(self, name, value):
        if name.startswith('_'):
            # Set private attributes directly on this proxy
            object.__setattr__(self, name, value)
        else:
            # Apply changes to the atom and record diff
            setattr(self._parent._atoms[self._index], name, value)
            diff = self._parent._capture_state_diff()
            # Send diff if it exists
            if diff and hasattr(self._parent, '_send_diff'):
                logger.info(f"posting diff from atom proxy: {diff}")
                self._parent._send_diff(diff)


class _PositionProxy:
    """Proxy for atom position that tracks changes"""
    
    def __init__(self, parent, index):
        self._parent = parent
        self._index = index
        self._array = parent._atoms.positions[index].copy()
    
    def _capture_and_send_diff(self):
        """Helper method to capture and send diffs"""
        diff = self._parent._capture_state_diff()
        # Send diff if it exists
        if diff and hasattr(self._parent, '_send_diff'):
            logger.info(f"posting diff from position proxy: {diff}")
            self._parent._send_diff(diff)

    def __array__(self, dtype=None):
        """Make this behave like a numpy array"""
        if dtype is not None:
            return np.array(self._array, dtype=dtype)
        return self._array
    
    def __getitem__(self, i):
        return self._parent._atoms.positions[self._index][i]
    
    def __setitem__(self, i, value):
        # Update both the local array and the actual atom position
        self._array[i] = value
        self._parent._atoms.positions[self._index][i] = value
        self._capture_and_send_diff()
    
    def __iadd__(self, other):
        # Handle += operation
        new_pos = self._parent._atoms.positions[self._index] + other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._capture_and_send_diff()
        return self
    
    def __isub__(self, other):
        # Handle -= operation
        new_pos = self._parent._atoms.positions[self._index] - other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._capture_and_send_diff()
        return self
    
    def __imul__(self, other):
        # Handle *= operation
        new_pos = self._parent._atoms.positions[self._index] * other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._capture_and_send_diff()
        return self
    
    def __itruediv__(self, other):
        # Handle /= operation
        new_pos = self._parent._atoms.positions[self._index] / other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._capture_and_send_diff()
        return self
    
    def __len__(self):
        return 3  # Positions are always 3D
    
    def copy(self):
        return self._parent._atoms.positions[self._index].copy()
    
    def tolist(self):
        return self._parent._atoms.positions[self._index].tolist()
