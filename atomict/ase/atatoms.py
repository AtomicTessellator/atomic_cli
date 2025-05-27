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
from atomict.api import patch, post


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
                logger.info(f"Capturing diff after {func.__name__}")
                self._capture_state_diff()
            
            return result
        return wrapper
    return decorator


class ATAtoms:
    """
    A lightweight wrapper around ASE Atoms that transparently tracks modifications
    and generates diffs for server synchronization.
    """
    
    def __init__(self, atoms: Atoms, server_url: Optional[str] = None, 
                 batch_size: int = 10, sync_interval: float = 10.0, project_id: Optional[str] = '',
                 batch_diffs: bool = False):
        """
        Initialize the ATAtoms wrapper.
        
        Parameters:
        -----------
        atoms: ASE Atoms object to wrap
        server_url: URL of the server to send diffs to (if None, will use AT_SERVER env var)
        batch_size: Number of diffs to accumulate before sending to server in a single request
        sync_interval: Maximum time in seconds between syncs to server (default: 10 seconds)
        project_id: ID of the project to associate with the atoms
        batch_diffs: Whether to batch diffs
        """
        if not isinstance(atoms, Atoms):
            raise TypeError(f"Expected ASE Atoms object, got {type(atoms).__name__}: {atoms}")

        self._atoms = atoms
        self._server_url = server_url or os.environ.get('AT_SERVER')
        self._project_id = project_id
        self._object_id = str(uuid.uuid4())
        self._batch_size = batch_size
        self._sync_interval = sync_interval
        self._diffs = []
        self._last_sync_time = time.time()
        self._seq_num = 0
        self._state_id = None
        self._run_id = None
        self._initialized = False
        self._batch_diffs = batch_diffs
        self._initial_state = self._get_current_state()
        self._previous_state = self._serialize_state(self._initial_state)
        self._structure_id = self._hash_state(self._previous_state)
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
        # If no server URL is set, skip the server initialization
        if not self._server_url:
            logger.info("No server URL provided. Skipping server initialization.")
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
                'api/atatoms-states/', 
                state_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            # Store the state ID for future reference
            if response and 'id' in response:
                self._state_id = response['id']
                print(f"Successfully initialized state on server with ID: {self._state_id}")
                
                # Initialize a run with this state
                if self._state_id:
                    self._initialize_run_on_server()
            else:
                print(f"Warning: Server response did not contain expected 'id' field: {response}")
            
        except Exception as e:
            print(f"Warning: Failed to initialize object on server: {e}")
    
    def _initialize_run_on_server(self):
        """Create an AtomicRun object on the server using the initial state"""
        # If no server URL is set, skip the run initialization
        if not self._server_url:
            logger.info("No server URL provided. Skipping run initialization.")
            return
            
        if not self._state_id:
            print("Warning: Cannot initialize run without a valid state ID")
            return
            
        try:
            # Prepare run data
            run_data = {
                'id': str(uuid.uuid4()),
                'name': f"ATAtomsRun-{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                'starting_state': self._state_id,
                'end_state': self._state_id,  # Initially same as starting state
                'project': self._project_id,
                'metadata': json.dumps({
                    'created_by': 'ATAtoms',
                    'initial_structure_id': self._structure_id,
                    'timestamp': datetime.datetime.now().isoformat()
                })
            }
            # Use the API module to create the run
            response = post(
                'api/atatoms-runs/',
                run_data,
            )
            
            # Store the run ID for future reference
            if response and 'id' in response:
                self._run_id = response['id']
                print(f"Successfully initialized run on server with ID: {self._run_id}")
            else:
                print(f"Warning: Server response did not contain expected 'id' field: {response}")
                
        except Exception as e:
            print(f"Warning: Failed to initialize run on server: {e}")
    
    def _capture_state_diff(self):
        """
        Capture differences between current and previous state
        """
        if not self._server_url:
            logger.info("No server configured, skipping creating diffs.")

        current_state = self._get_current_state()
        serialized_current = self._serialize_state(current_state)
        
        # For initial state, just initialize on server without creating a diff
        if not self._initialized:
            logger.info("Initializing on server (first state)")
            self._initialize_on_server()
            self._initialized = True
            return None
        
        diff = DeepDiff(self._previous_state, serialized_current, verbose_level=1)
        
        if diff:
            timestamp = datetime.datetime.now().isoformat()
            logger.info(f"State change detected, creating diff with seq_num={self._seq_num}")
            diff_item = {
                'timestamp': timestamp,
                'sequence': self._seq_num,
                'diff': diff
            }
            
            self._previous_state = serialized_current
            self._seq_num += 1
            logger.info(f"Incremented internal seq_num to {self._seq_num}")
            
            if self._batch_diffs:
                # Add to batch if batching is enabled
                logger.info(f"Adding diff to batch queue (queue size now: {len(self._diffs)+1})")
                self._diffs.append(diff_item)
                
                # Check if we should sync based on batch size or time interval
                current_time = time.time()
                if (len(self._diffs) >= self._batch_size or 
                    current_time - self._last_sync_time >= self._sync_interval):
                    logger.info(f"Batch threshold reached, syncing {len(self._diffs)} diffs")
                    self._sync_diffs()
            else:
                self._send_diff(diff)
        else:
            logger.info("No state change detected")
    
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

    def _recursive_serialize(self, obj):
        """Recursively convert non-serializable objects to serializable forms."""
        import numpy as np
        # Handle numpy arrays
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        # Handle numpy scalars
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj) if isinstance(obj, np.floating) else int(obj)
        # Handle dicts
        elif isinstance(obj, dict):
            return {k: self._recursive_serialize(v) for k, v in obj.items()}
        # Handle lists/tuples
        elif isinstance(obj, (list, tuple)):
            return [self._recursive_serialize(v) for v in obj]
        # Handle objects with .todict() (e.g., Spacegroup)
        elif hasattr(obj, 'todict') and callable(getattr(obj, 'todict')):
            return self._recursive_serialize(obj.todict())
        # Handle objects with __dict__ (as a fallback, but avoid recursion on builtins)
        elif hasattr(obj, '__dict__') and not isinstance(obj, type):
            return {k: self._recursive_serialize(v) for k, v in obj.__dict__.items() if not k.startswith('__')}
        # Handle basic types
        elif isinstance(obj, (str, int, float, bool)) or obj is None:
            return obj
        # Fallback: string representation
        else:
            return str(obj)

    def _serialize_state(self, state: Dict[str, Any]) -> Dict[str, Any]:
        """Convert all non-serializable objects in the state to serializable forms recursively."""
        return self._recursive_serialize(state)
    
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
            # Update our internal atoms reference while preserving tracking state
            self._atoms = result
            # Capture the state change
            self._capture_state_diff()
            return self
            
    # Explicitly track rattle which modifies positions directly
    @track_state_changes()
    def rattle(self, *args, **kwargs):
        result = self._atoms.rattle(*args, **kwargs)
        return self if result is self._atoms else result
    
    # Handle individual atom access with change tracking
    def __getitem__(self, index):
        if isinstance(index, slice):
            # Get the sliced atoms
            new_atoms = self._atoms[index]
            
            # Instead of creating a new ATAtoms instance, create a copy of self
            # and update its _atoms attribute to the sliced atoms
            if new_atoms is self._atoms:  # If the slice didn't change anything
                return self
            else:
                # Update internal atoms but preserve all tracking state
                self._atoms = new_atoms
                # Record the state change
                self._capture_state_diff()
                return self
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
    
    # Modify __getattr__ to handle methods that return None but modify state
    def __getattr__(self, name):
        attr = getattr(self._atoms, name)
        
        if callable(attr):
            def wrapped_method(*args, **kwargs):
                result = attr(*args, **kwargs)
                if result is not None:
                    logger.info(f"Result: {result}")
                # Handle methods that return the atoms object
                if result is self._atoms:
                    self._capture_state_diff()
                    return self
                elif isinstance(result, Atoms):
                    return ATAtoms(result)
                
                self._capture_state_diff()
                return self if result is None else result

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
        # Skip if no server URL is set
        if not self._server_url:
            logger.info("No server URL provided. Skipping diff send.")
            return None
            
        # Don't attempt to send if we don't have a state_id yet
        if not self._state_id:
            logger.warning("Cannot send diff without state_id. Attempting to initialize on server first.")
            self._initialize_on_server()
            if not self._state_id:
                return None
        
        try:
            # Get current state to include in the request
            # current_state = self._get_current_state()
            # serialized_current = self._serialize_state(current_state)
            
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
                'state': self._state_id,  # Use state_id as the reference to the state
                'run': self._run_id  # Include the run ID if available
            }
            
            logger.info(f"Sending diff to server with sequence={self._seq_num - 1}")
            
            # Send the diff to the server
            response = post(
                'api/atatoms-diffs/', 
                diff_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            if not response:
                logger.warning("Empty response when sending diff to server")
            
            return response
        except Exception as e:
            logger.error(f"Failed to send diff to server: {e}")
            return None

    def save_current_state(self) -> Dict[str, Any]:
        """
        Get current state, serialize it, hash it, and send it to server.
        
        Returns:
        --------
        Dict with server response data
        """
        # Skip if no server URL is set
        if not self._server_url:
            logger.info("No server URL provided. Skipping state save.")
            return {"info": "No server URL provided. State was not saved."}
            
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
                'api/atatoms-states/', 
                state_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            if response and 'id' in response:
                logger.info(f"Successfully saved state on server with ID: {response['id']}")
                
                # Update the end state of the run if available
                if self._run_id:
                    self.update_run_end_state(response['id'])
                    
                return response
            else:
                logger.warning(f"Server response did not contain expected 'id' field: {response}")
                return response
                
        except Exception as e:
            logger.error(f"Failed to save state on server: {e}")
            return {"error": str(e)}

    def update_run_end_state(self, state_id):
        """Update the end state of the current run"""
        if not self._server_url:
            logger.info("No server URL provided. Skipping run end state update.")
            return
            
        if not self._run_id:
            logger.warning("Cannot update run end state: no run ID available")
            return
            
        try:
            # Prepare update data
            update_data = {
                'end_state': state_id,
                'updated_at': datetime.datetime.now().isoformat()
            }
            
            # Use the API module to update the run
            response = patch(
                f'api/atatoms-runs/{self._run_id}/',
                update_data,
            )
            
            if response and 'id' in response:
                logger.info(f"Successfully updated run end state to {state_id}")
            else:
                logger.warning(f"Failed to update run end state: {response}")
                
        except Exception as e:
            logger.error(f"Failed to update run end state: {e}")

    @classmethod
    def from_known_state(cls, state_id: str) -> 'ATAtoms':
        """Retrieves a state from the server and initializes an ATAtoms object with it."""
        raise NotImplementedError

    def track_changes(self):
        """
        DEPRECATED

        Simple callback for ASE optimizers to track changes.
        Calculates diff between previous and current state and sends if changes exist.
        
        Usage:
            opt = BFGS(atoms)
            opt.attach(atoms.track_changes, interval=1)
            opt.run(fmax=0.02)
        """
        logger.info("track_changes callback invoked")
        self._capture_state_diff()
        # this may work without attaching now
        return True

    def _sync_diffs(self):
        """Send accumulated diffs to the server in a single request and clear the queue"""
        # Skip if no server URL is set
        if not self._server_url:
            logger.info("No server URL provided. Skipping diffs sync and pruning.")
            self._diffs = []
            return
            
        if not self._diffs:
            return
        
        if not self._state_id:
            self._initialize_on_server()
            if not self._state_id:
                return
            
        try:
            # Prepare batch of diffs
            batch_data = {
                'atomic_state_id': self._state_id,
                'structure_id': self._structure_id,
                'run': self._run_id,
                'diffs': []
            }
            
            # Add each diff to the batch
            for diff_item in self._diffs:
                serialized_diff = self._serialize_diff(diff_item['diff'])
                
                batch_data['diffs'].append({
                    'sequence': diff_item['sequence'],
                    'timestamp': diff_item['timestamp'],
                    'data': serialized_diff
                })
            
            # Send the batch in a single request
            response = post(
                'api/atatoms-diffs/batch/',
                batch_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            logger.info(f"Sent batch of {len(self._diffs)} diffs to server")
            self._diffs = []
            self._last_sync_time = time.time()
            
        except Exception as e:
            logger.error(f"Failed to sync diff batch: {e}")


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
        logger.info(f"Setting position for atom {self._index}")
        self._parent._atoms.positions[self._index] = value
        self._parent._capture_state_diff()
    
    def __getattr__(self, name):
        # Forward all other attributes to the actual atom
        return getattr(self._parent._atoms[self._index], name)
    
    def __setattr__(self, name, value):
        if name.startswith('_'):
            # Set private attributes directly on this proxy
            object.__setattr__(self, name, value)
        else:
            # Apply changes to the atom and record diff
            logger.info(f"Setting attribute {name} for atom {self._index}")
            setattr(self._parent._atoms[self._index], name, value)
            self._parent._capture_state_diff()


class _PositionProxy:
    """Proxy for atom position that tracks changes"""
    
    def __init__(self, parent, index):
        self._parent = parent
        self._index = index
        self._array = parent._atoms.positions[index].copy()
    
    def _capture_and_send_diff(self):
        """Helper method to capture and send diffs"""
        logger.info(f"Capturing diff from position proxy for atom {self._index}")
        self._parent._capture_state_diff()

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
