import uuid
import time
import numpy as np
import json
import datetime
import requests
import os
import sys
import functools
from typing import Optional, Dict, Any, List, Union, Tuple
from ase.atoms import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT

from deepdiff import DeepDiff
# t1 = {1:1, 2:2, 3:3}

# t2 = {1:1, 2:4, 3:3}

# pprint(DeepDiff(t1, t2, verbose_level=0), indent=2)
# {'values_changed': {'root[2]': {'new_value': 4, 'old_value': 2}}}

def track_state_changes():
    """
    Decorator to track state changes after method execution.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            # Call the original method
            result = func(self, *args, **kwargs)
            
            # Track state changes after the method call
            if hasattr(self, '_capture_state_diff'):
                self._capture_state_diff()
            
            return result
        return wrapper
    return decorator

class ATAtoms:
    """
    A lightweight wrapper around ASE Atoms that transparently tracks modifications
    and generates deltas for server synchronization.
    """
    
    def __init__(self, atoms: Atoms, server_url: Optional[str] = None, 
                 batch_size: int = 10, sync_interval: float = 60.0):
        """
        Initialize the ATAtoms wrapper.
        
        Parameters:
        -----------
        atoms: ASE Atoms object to wrap
        server_url: URL of the server to send deltas to (if None, will use ATOMICT_SERVER_URL env var)
        batch_size: Number of deltas to accumulate before sending to server
        sync_interval: Maximum time in seconds between syncs to server
        """
        # Initialize core properties
        self._atoms = atoms
        self._object_id = str(uuid.uuid4())
        self._server_url = server_url or os.environ.get('ATOMICT_SERVER_URL')
        self._batch_size = batch_size
        self._sync_interval = sync_interval
        self._deltas = []
        self._diffs = []  # Track the state diffs
        self._last_sync_time = time.time()
        self._execution_id = str(uuid.uuid4())
        self._token = os.environ.get('ATOMICT_TOKEN', "")
        self._seq_num = 0  # Initialize sequence number counter
        
        # Initialize state tracking
        self._initial_state = self._get_current_state()  # Save the initial state
        self._previous_state = self._serialize_state(self._initial_state)  # Serialize for comparison
        
        # Record the initial state
        self._capture_state_diff()
    
    def _capture_state_diff(self):
        """
        Capture differences between current and previous state and append to diffs array
        """
        current_state = self._get_current_state()
        serialized_current = self._serialize_state(current_state)
        
        # For initial state, just store the full state
        if not self._diffs:
            self._diffs.append({
                'timestamp': datetime.datetime.now().isoformat(),
                'sequence': self._seq_num,
                'state': serialized_current
            })
            
            # Also record as a delta for server sync
            self._deltas.append({
                'object_id': self._object_id,
                'execution_id': self._execution_id,
                'timestamp': datetime.datetime.now().isoformat(),
                'sequence': self._seq_num,
                'state': serialized_current
            })
            
            self._previous_state = serialized_current
            self._seq_num += 1  # Increment sequence number
            return
        
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
            
            # Also record as a delta for server sync
            self._deltas.append({
                'object_id': self._object_id,
                'execution_id': self._execution_id,
                'timestamp': timestamp,
                'sequence': self._seq_num,
                'diff': diff
            })
            
            # Update previous state for next comparison
            self._previous_state = serialized_current
            self._seq_num += 1  # Increment sequence number
            
            # Check if we should sync now
            if len(self._deltas) >= self._batch_size or \
               (time.time() - self._last_sync_time) >= self._sync_interval:
                self._sync_deltas()
    
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
            if atoms.positions is not None:
                state['positions'] = atoms.positions.copy()
                
            # Handle cell
            if hasattr(atoms.cell, 'array'):
                state['cell'] = atoms.cell.array.copy()
            else:
                state['cell'] = atoms.cell.copy() if hasattr(atoms.cell, 'copy') else atoms.cell
                
            # Handle pbc
            state['pbc'] = atoms.pbc.copy() if hasattr(atoms.pbc, 'copy') else atoms.pbc
            
            # Handle atomic numbers and symbols
            state['numbers'] = atoms.numbers.copy()
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
        
        # Add any arrays that might have been set
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
    
    def _sync_deltas(self):
        """Send accumulated deltas to the server"""
        if not self._server_url or not self._deltas:
            return
        
        try:
            headers = {
                'Content-Type': 'application/json',
                'Accept': 'application/json',
                'User-Agent': 'ATAtoms/1.0'
            }
            
            if self._token:
                headers['Authorization'] = f"Token {self._token}"
            
            for delta in self._deltas:
                response = requests.post(
                    f"{self._server_url}/api/atatoms-delta/",
                    data=json.dumps(delta, ensure_ascii=False),
                    headers=headers
                )
                
                if response.status_code >= 400:
                    print(f"Warning: Server returned error {response.status_code}: {response.text}")
            
            self._deltas = []
            self._last_sync_time = time.time()
            
        except Exception as e:
            print(f"Warning: Failed to sync deltas with server: {e}")
    
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
        """Ensure any remaining deltas are synced before garbage collection"""
        if hasattr(self, '_deltas') and self._deltas:
            try:
                self._sync_deltas()
            except:
                pass


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
        self._parent._capture_state_diff()
    
    def __getattr__(self, name):
        # Forward all other attributes to the actual atom
        return getattr(self._parent._atoms[self._index], name)
    
    def __setattr__(self, name, value):
        if name.startswith('_'):
            # Set private attributes directly on this proxy
            object.__setattr__(self, name, value)
        else:
            # Apply changes to the atom and record delta
            setattr(self._parent._atoms[self._index], name, value)
            self._parent._capture_state_diff()


class _PositionProxy:
    """Proxy for atom position that tracks changes"""
    
    def __init__(self, parent, index):
        self._parent = parent
        self._index = index
        self._array = parent._atoms.positions[index].copy()
    
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
        self._parent._capture_state_diff()
    
    def __iadd__(self, other):
        # Handle += operation
        new_pos = self._parent._atoms.positions[self._index] + other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._parent._capture_state_diff()
        return self
    
    def __isub__(self, other):
        # Handle -= operation
        new_pos = self._parent._atoms.positions[self._index] - other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._parent._capture_state_diff()
        return self
    
    def __imul__(self, other):
        # Handle *= operation
        new_pos = self._parent._atoms.positions[self._index] * other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._parent._capture_state_diff()
        return self
    
    def __itruediv__(self, other):
        # Handle /= operation
        new_pos = self._parent._atoms.positions[self._index] / other
        self._parent._atoms.positions[self._index] = new_pos
        self._array = new_pos.copy()
        self._parent._capture_state_diff()
        return self
    
    def __len__(self):
        return 3  # Positions are always 3D
    
    def copy(self):
        return self._parent._atoms.positions[self._index].copy()
    
    def tolist(self):
        return self._parent._atoms.positions[self._index].tolist()
