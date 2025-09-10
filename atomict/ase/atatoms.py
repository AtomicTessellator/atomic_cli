import uuid
import time
import numpy as np
import json
import datetime
import os
import sys
import hashlib
import threading
import queue
from typing import Optional, Dict, Any
import logging

from ase.atoms import Atoms
from deepdiff import DeepDiff

from atomict.api import patch, post


logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.ERROR)


class ATAtomsError(Exception):
    pass


class ATAtomsServerError(ATAtomsError):
    """Raised when server operations fail"""
    pass


class ATAtomsInitializationError(ATAtomsServerError):
    """Raised when server initialization fails"""
    pass


class QueueWorker:
    """Handles network operations in a separate thread"""
    def __init__(self):
        self._queue = queue.Queue(maxsize=100_000)
        self._thread = threading.Thread(target=self._process_queue, daemon=True)
        self._thread.start()
        
    def _process_queue(self):
        """Process network operations from the queue"""
        while True:
            try:
                operation = self._queue.get()
                if operation is None:
                    break
                    
                func, args, kwargs = operation
                try:
                    func(*args, **kwargs)
                except Exception as e:
                    logger.error(f"Network operation failed: {e}")
                finally:
                    self._queue.task_done()
            except Exception as e:
                logger.error(f"Error processing network queue: {e}")
                
    def enqueue(self, func, *args, **kwargs):
        """Add a request to the queue without blocking"""
        try:
            self._queue.put((func, args, kwargs), block=False)
        except queue.Full:
            logger.warning("Network operation queue is full, operation dropped")


class ATAtoms:
    """
    A wrapper around ASE Atoms that transparently tracks modifications
    and tracks diffs between Atoms states. An ATAtoms run is a collection of diffs.
    """
    
    def __init__(self, 
                 atoms: Atoms, 
                 project_id: Optional[str] = '',
                 simulation_id: Optional[str] = '', 
                 batch_size: int = 20, 
                 sync_interval: float = 10.0, 
                 batch_diffs: bool = False, 
                 is_result: bool = False):
        """
        atoms: ASE Atoms object to wrap
        project_id: ID of the project to associate with the atoms (env var AT_PROJECT_ID)
        simulation_id: ID of simulation the object is being used in (env var AT_SIMULATION_ID)
        batch_size: Number of diffs to accumulate before sending to server in a single request
        sync_interval: Maximum time in seconds between syncs to server
        batch_diffs: Whether to batch diffs
        is_result: A (metadata) flag for indicating whether the object is the result of a simulation
        """
        if not isinstance(atoms, Atoms):
            raise TypeError(f"Expected ASE Atoms object, got {type(atoms).__name__}: {atoms}")

        self._atoms = atoms
        self._project_id = project_id or os.environ.get('AT_PROJECT_ID')
        self._simulation_id = simulation_id or os.environ.get('AT_SIMULATION_ID')
        if not self._project_id or not self._simulation_id:
            raise ATAtomsInitializationError("ATAtoms requires project_id and simulation_id to be passed or set via environment variables (AT_PROJECT_ID, AT_SIMULATION_ID)")

        self._object_id = str(uuid.uuid4())
        self._batch_size = batch_size
        self._sync_interval = sync_interval
        self._worker = self._create_worker()
        self._diffs = []
        self._last_sync_time = time.time()
        self._seq_num = 0
        self.atomic_state_id = None
        self._run_id = None
        self._initialized = False
        self._batch_diffs = batch_diffs
        self._is_result = is_result
        self._initial_state = self._get_current_state()
        self._previous_state = self._serialize_state(self._initial_state)
        self._structure_id = self._hash_state(self._previous_state)
        self._capture_state_diff()
    
    def _create_worker(self) -> QueueWorker:
        """Create and return a new QueueWorker instance.
        
        This method is separated out to allow for easier testing and mocking.
        """
        return QueueWorker()
    
    def _hash_state(self, state_data):
        """Generate a SHA-256 hash of the state data"""
        state_json = json.dumps(state_data, sort_keys=True)
        _hash = hashlib.sha256(state_json.encode('utf-8')).hexdigest()
        logger.info(f"Created structure id: {_hash}")
        return _hash
    
    def _initialize_on_server(self):
        """Send the initial state to create an object on the server"""
        # If no server URL is set, skip the server initialization and run offline
        
        try:
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
            response = post(
                'api/atatoms-states/', 
                state_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            if response and 'id' in response:
                self.atomic_state_id = response['id']
                logger.info(f"Successfully initialized state on server with ID: {self.atomic_state_id}")
                if self.atomic_state_id:
                    self._initialize_run_on_server()
            else:
                logger.warning(f"Server response did not contain expected 'id' field: {response}")
            
        except Exception as e:
            logger.error(f"Failed to queue initialization on server: {e}")
    
    def _initialize_run_on_server(self):
        """Create an AtomicRun object on the server using the initial state"""
        # If no server URL is set, skip the run initialization and run offline
            
        if not self.atomic_state_id:
            logger.warning("Cannot initialize run without a valid state ID")
            return
            
        try:
            run_data = {
                'id': str(uuid.uuid4()),
                'name': f"ATAtomsRun-{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                'starting_state': self.atomic_state_id,
                'end_state': self.atomic_state_id,
                'project': self._project_id,
                'metadata': json.dumps({
                    'created_by': 'ATAtoms',
                    'initial_structure_id': self._structure_id,
                    'timestamp': datetime.datetime.now().isoformat(),
                    'is_result': self._is_result,
                    'simulation_id': self._simulation_id
                })
            }
            response = post(
                'api/atatoms-runs/',
                run_data,
            )
            if response and 'id' in response:
                self._run_id = response['id']
                logger.info(f"Successfully initialized run on server with ID: {self._run_id}")
            else:
                logger.warning(f"Server response did not contain expected 'id' field: {response}")
                
        except Exception as e:
            logger.error(f"Failed to queue run initialization on server: {e}")
    
    def _capture_state_diff(self):
        """
        Capture differences between current and previous state
        """

        current_state = self._get_current_state()
        serialized_current = self._serialize_state(current_state)
        new_atomic_state_id = self._hash_state(serialized_current)
        if self._structure_id == new_atomic_state_id:
            logger.info("New state hash matches current")
            return

        if not self._initialized:
            logger.info("Initializing on server (first state)")
            self._initialize_on_server()
            self._initialized = True
            return None
        
        diff = DeepDiff(self._previous_state, serialized_current, verbose_level=1)
        
        if diff:
            timestamp = datetime.datetime.now().isoformat()
            current_seq = self._seq_num
            logger.info(f"State change detected, creating diff with seq_num={current_seq}")
            diff_item = {
                'timestamp': timestamp,
                'sequence': current_seq,
                'diff': diff
            }
            
            self._previous_state = serialized_current
            self._structure_id = self._hash_state(self._previous_state)
            
            if self._batch_diffs:
                logger.info(f"Adding diff to batch queue (queue size now: {len(self._diffs)+1})")
                self._diffs.append(diff_item)
                current_time = time.time()
                if (len(self._diffs) >= self._batch_size or 
                    current_time - self._last_sync_time >= self._sync_interval):
                    logger.info(f"Batch threshold reached, syncing {len(self._diffs)} diffs")
                    self._sync_diffs()
            else:
                self._send_diff(diff)

            self._seq_num += 1
            logger.info(f"Incremented internal seq_num to {self._seq_num}")
        else:
            logger.info("No state change detected")
    
    @property
    def calc(self):
        return self._atoms.calc
    
    @calc.setter
    def calc(self, calculator):
        self._atoms.calc = calculator
        self._capture_state_diff()
    
    def _get_current_state(self) -> Dict[str, Any]:
        """Get the complete state of the atoms object"""
        state = self._atoms.todict()

        if 'symbols' not in state:
            state['symbols'] = self._atoms.get_chemical_symbols()

        if hasattr(self._atoms, 'get_potential_energy'):
            try:
                state['energy'] = self._atoms.get_potential_energy()
            except:
                pass
                
        if hasattr(self._atoms, 'get_forces'):
            try:
                state['forces'] = self._atoms.get_forces()
            except:
                pass
                
        if hasattr(self._atoms, 'get_stress'):
            try:
                state['stress'] = self._atoms.get_stress()
            except:
                pass

        if hasattr(self._atoms, 'info') and self._atoms.info:
            state['info'] = self._atoms.info.copy()

        for name in self._atoms.arrays:
            if name not in ['positions', 'numbers']:
                state[name] = self._atoms.arrays[name].copy()
                
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

        serialized_diff = {}        

        for diff_type, values in diff.items():
            if hasattr(values, 'items'):
                serialized_diff[diff_type] = dict(values)
            elif isinstance(values, (list, tuple)):
                serialized_diff[diff_type] = list(values)
            elif hasattr(values, '__iter__') and not isinstance(values, str):
                serialized_diff[diff_type] = list(values)
            else:
                serialized_diff[diff_type] = values
                
        return serialized_diff
    
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
                self._capture_state_diff()
                return self
        else:
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

    def __getattr__(self, name):
        # First check if the attribute exists directly on the class or instance
        if name in dir(self.__class__) or name in self.__dict__:
            return object.__getattribute__(self, name)
            
        # If not found on self, try to get it from the wrapped atoms object
        try:
            attr = getattr(self._atoms, name)
        except AttributeError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
        
        if callable(attr):
            def wrapped_method(*args, **kwargs):
                result = attr(*args, **kwargs)
                if result is self._atoms:
                    self._capture_state_diff()
                    return self
                elif isinstance(result, Atoms):
                    self._atoms = result
                    self._capture_state_diff()
                    return self
            
                self._capture_state_diff()
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
                self._worker._queue.join(timeout=5.0)
            except:
                pass

    def __enter__(self):
        """Context manager entry point"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point - syncs diffs and cleans up worker thread"""
        try:
            if hasattr(self, '_diffs') and self._diffs:
                self._sync_diffs()

            if hasattr(self, '_worker'):
                self._worker._queue.join()

            if hasattr(self, '_worker'):
                self._worker._queue.put(None)
                self._worker._thread.join(timeout=1.0)
                
        except Exception as e:
            logger.error(f"Error during ATAtoms cleanup: {e}")
            return False
            
        return False

    def _send_diff(self, diff):
        """
        Send a diff to the server at /api/atatoms-diffs
        """

        if not self.atomic_state_id:
            logger.warning("Cannot send diff without state_id. Attempting to initialize on server first.")
            self._initialize_on_server()
            if not self.atomic_state_id:
                raise ATAtomsInitializationError("Failed to initialize ATAtoms state on server after retry attempt")

        try:
            serialized_diff = self._serialize_diff(diff)
            current_seq = self._seq_num
            diff_data = {
                'atoms_id': self._object_id,
                'sequence': current_seq,
                'timestamp': datetime.datetime.now().isoformat(),
                'data': serialized_diff,
                'run': self._run_id
            }
            
            logger.info(f"Queueing diff to server with sequence={current_seq}")
            
            def _do_send_diff(diff_data):
                response = post(
                    'api/atatoms-diffs/', 
                    diff_data,
                    extra_headers={'Content-Type': 'application/json'}
                )
                
                if not response:
                    logger.warning("Empty response when sending diff to server")

            self._worker.enqueue(_do_send_diff, diff_data)
            
        except Exception as e:
            logger.error(f"Failed to queue diff send to server: {e}")
            return None

    def save_current_state(self) -> Dict[str, Any]:
        """
        Get current state, serialize it, hash it, and send it to server.
        
        Returns:
        --------
        Dict with server response data
        """

        current_state = self._get_current_state()
        serialized_state = self._serialize_state(current_state)
        structure_id = self._hash_state(serialized_state)
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
            logger.info(f"Saving current state with structure_id: {structure_id}")
            response = post(
                'api/atatoms-states/', 
                state_data,
                extra_headers={'Content-Type': 'application/json'}
            )
            
            if response and 'id' in response:
                logger.info(f"Successfully saved state on server with ID: {response['id']}")
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
            
        if not self._run_id:
            logger.warning("Cannot update run end state: no run ID available")
            return
            
        try:
            update_data = {
                'end_state': state_id,
                'updated_at': datetime.datetime.now().isoformat()
            }
            
            def _do_update_run():
                response = patch(
                    f'api/atatoms-runs/{self._run_id}/',
                    update_data,
                )
                
                if response and 'id' in response:
                    logger.info(f"Successfully updated run end state to {state_id}")
                else:
                    logger.warning(f"Failed to update run end state: {response}")

            self._worker.enqueue(_do_update_run)
                
        except Exception as e:
            logger.error(f"Failed to queue run end state update: {e}")

    @classmethod
    def from_known_state(cls, state_id: str) -> 'ATAtoms':
        """Retrieves a state from the server and initializes an ATAtoms object with it."""
        raise NotImplementedError

    def _sync_diffs(self):
        """Send accumulated diffs to the server in a single request and clear the queue"""
        # offline mode
            
        if not self._diffs:
            return
        
        if not self.atomic_state_id:
            self._initialize_on_server()
            if not self.atomic_state_id:
                return
            
        try:
            batch_data = {
                'run': self._run_id,
                'diffs': []
            }
            for diff_item in self._diffs:
                serialized_diff = self._serialize_diff(diff_item['diff'])
                
                batch_data['diffs'].append({
                    'sequence': diff_item['sequence'],
                    'timestamp': diff_item['timestamp'],
                    'data': serialized_diff
                })

            self._diffs = []
            self._last_sync_time = time.time()
            
            def _do_sync_diffs(batch_data):
                _ = post(
                    'api/atatoms-diffs/batch/',
                    batch_data,
                    extra_headers={'Content-Type': 'application/json'}
                )
                logger.info(f"Sent batch of {len(batch_data['diffs'])} diffs to server")

            self._worker.enqueue(_do_sync_diffs, batch_data)
            
        except Exception as e:
            logger.error(f"Failed to queue diff batch sync: {e}")


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
