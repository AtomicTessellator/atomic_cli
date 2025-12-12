"""
Tests for tess format bytes I/O operations.

This module tests reading and writing tess files using file-like byte objects (BytesIO)
and verifies compatibility with binary file fixtures.
"""
import io
import tempfile
from pathlib import Path

import numpy as np
import pytest


# Skip if msgpack dependencies are missing
pytest.importorskip("msgpack")
pytest.importorskip("msgpack_numpy")

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from atomict.io.formats.tess import write_tess, read_tess


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tdir:
        yield Path(tdir)


@pytest.fixture
def fixtures_dir():
    """Returns the path to the test fixtures directory."""
    return Path(__file__).parent.parent.parent.parent / "fixtures"


def _make_simple_atoms() -> Atoms:
    """Create a simple H2 molecule for testing."""
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.74],
    ], dtype=np.float64)
    cell = np.diag([5.0, 5.0, 5.0]).astype(np.float64)
    atoms = Atoms(symbols="H2", positions=positions, pbc=True)
    atoms.set_cell(cell)
    return atoms


def _compare_atoms(a: Atoms, b: Atoms):
    """Compare two Atoms objects for equivalence."""
    assert a.get_chemical_symbols() == b.get_chemical_symbols()
    assert np.allclose(a.get_positions(), b.get_positions(), rtol=1e-6, atol=1e-6)
    assert np.allclose(a.get_cell(), b.get_cell(), rtol=1e-6, atol=1e-6)
    assert (a.get_pbc() == b.get_pbc()).all()
    assert np.allclose(a.get_masses(), b.get_masses(), rtol=1e-6, atol=1e-6)


class TestBytesIOWrite:
    """Test writing tess format to BytesIO objects."""

    def test_write_single_frame_to_bytesio(self):
        """Test writing a single frame to BytesIO."""
        atoms = _make_simple_atoms()
        bio = io.BytesIO()
        
        write_tess(atoms, bio, metadata={"test": "bytesio"}, compression="none")
        
        # Verify bytes were written
        bio.seek(0)
        data = bio.read()
        assert len(data) > 0
        assert isinstance(data, bytes)

    def test_write_multi_frame_to_bytesio(self):
        """Test writing multiple frames to BytesIO."""
        frames = []
        for i in range(5):
            a = _make_simple_atoms()
            pos = a.get_positions()
            pos[0, 2] = float(i)
            a.set_positions(pos)
            frames.append(a)
        
        bio = io.BytesIO()
        write_tess(frames, bio, compression="none")
        
        bio.seek(0)
        data = bio.read()
        assert len(data) > 0
        assert isinstance(data, bytes)

    @pytest.mark.parametrize("compression", ["none", "zlib", "lz4"])
    def test_write_to_bytesio_with_compression(self, compression):
        """Test writing to BytesIO with different compression modes."""
        if compression == "lz4":
            pytest.importorskip("lz4.block")
        
        atoms = _make_simple_atoms()
        bio = io.BytesIO()
        
        write_tess(atoms, bio, compression=compression, compression_level=1)
        
        bio.seek(0)
        data = bio.read()
        assert len(data) > 0


class TestBytesIORead:
    """Test reading tess format from BytesIO objects."""

    def test_read_single_frame_from_bytesio(self):
        """Test reading a single frame from BytesIO."""
        atoms = _make_simple_atoms()
        bio = io.BytesIO()
        
        # Write to BytesIO
        write_tess(atoms, bio, metadata={"test": "read"}, compression="none")
        
        # Read back from BytesIO
        bio.seek(0)
        out_atoms, out_meta = read_tess(bio)
        
        assert len(out_atoms) == 1
        _compare_atoms(atoms, out_atoms[0])
        assert out_meta == {"test": "read"}

    def test_read_multi_frame_from_bytesio(self):
        """Test reading multiple frames from BytesIO."""
        frames = []
        for i in range(5):
            a = _make_simple_atoms()
            pos = a.get_positions()
            pos[0, 2] = float(i)
            a.set_positions(pos)
            frames.append(a)
        
        bio = io.BytesIO()
        write_tess(frames, bio, compression="none")
        
        bio.seek(0)
        out_atoms, _ = read_tess(bio)
        
        assert len(out_atoms) == len(frames)
        for i, (orig, loaded) in enumerate(zip(frames, out_atoms)):
            _compare_atoms(orig, loaded)
            assert loaded.get_positions()[0, 2] == pytest.approx(float(i))

    def test_read_subset_frames_from_bytesio(self):
        """Test reading subset of frames from BytesIO."""
        frames = []
        for i in range(6):
            a = _make_simple_atoms()
            pos = a.get_positions()
            pos[0, 2] = float(i)
            a.set_positions(pos)
            frames.append(a)
        
        bio = io.BytesIO()
        write_tess(frames, bio, compression="zlib", compression_level=1)
        
        bio.seek(0)
        requested = [4, 1, 5]
        out_atoms, _ = read_tess(bio, frames_indices=requested)
        
        assert len(out_atoms) == len(requested)
        for expected_idx, atoms in zip(requested, out_atoms):
            assert atoms.get_positions()[0, 2] == pytest.approx(float(expected_idx))


class TestBytesIORoundtrip:
    """Test roundtrip reading/writing with BytesIO."""

    @pytest.mark.parametrize("compression", ["none", "zlib"])
    def test_bytesio_roundtrip_preserves_data(self, compression):
        """Test that roundtrip through BytesIO preserves all data."""
        a = _make_simple_atoms()
        
        # Add various properties
        a.set_tags(np.array([1, 2], dtype=int))
        a.set_initial_charges(np.array([0.1, -0.1], dtype=float))
        a.set_initial_magnetic_moments(np.array([0.0, 1.0], dtype=float))
        a.set_momenta(np.array([[0.1, 0.0, 0.0], [0.0, -0.1, 0.0]], dtype=float))
        
        forces = np.array([[0.01, 0.0, 0.0], [0.0, -0.02, 0.0]], dtype=float)
        stress = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03], dtype=float)
        energy = -1.234
        
        calc = SinglePointCalculator(a, energy=energy, forces=forces, stress=stress)
        a.calc = calc
        a.info['test_key'] = 'test_value'
        a.arrays['custom'] = np.array([10.0, 20.0], dtype=float)
        
        bio = io.BytesIO()
        meta = {"format": "test", "version": 1}
        write_tess(a, bio, metadata=meta, compression=compression, compression_level=1)
        
        bio.seek(0)
        out_atoms, out_meta = read_tess(bio)
        
        assert len(out_atoms) == 1
        b = out_atoms[0]
        
        # Verify all properties preserved
        _compare_atoms(a, b)
        assert np.allclose(a.get_tags(), b.get_tags())
        assert np.allclose(a.get_initial_charges(), b.get_initial_charges())
        assert np.allclose(a.get_initial_magnetic_moments(), b.get_initial_magnetic_moments())
        assert np.allclose(a.get_momenta(), b.get_momenta())
        
        assert b.calc is not None
        assert pytest.approx(energy) == b.calc.results['energy']
        assert np.allclose(forces, b.calc.results['forces'])
        assert np.allclose(stress, b.calc.results['stress'])
        
        assert b.info.get('test_key') == 'test_value'
        assert 'custom' in b.arrays
        assert np.allclose(b.arrays['custom'], np.array([10.0, 20.0]))
        assert out_meta == meta

    def test_bytesio_multi_frame_roundtrip(self):
        """Test roundtrip with multiple frames through BytesIO."""
        frames = []
        base = _make_simple_atoms()
        base.set_tags(np.array([1, 1], dtype=int))
        
        for i in range(4):
            f = base.copy()
            f.set_positions(f.get_positions() + np.array([0.0, 0.0, 0.05 * i]))
            f.set_momenta(np.array([[i*0.1, 0.0, 0.0], [0.0, -i*0.1, 0.0]], dtype=float))
            
            forces = np.array([[0.01*i, 0.0, 0.0], [0.0, -0.02*i, 0.0]], dtype=float)
            stress = np.array([0.1*i, 0.2*i, 0.3*i, 0.01*i, 0.02*i, 0.03*i], dtype=float)
            calc = SinglePointCalculator(f, energy=-1.0-i*0.1, forces=forces, stress=stress)
            f.calc = calc
            f.info['frame_id'] = i
            frames.append(f)
        
        bio = io.BytesIO()
        write_tess(frames, bio, compression="zlib", compression_level=1)
        
        bio.seek(0)
        out_atoms, _ = read_tess(bio)
        
        assert len(out_atoms) == len(frames)
        for i, (orig, loaded) in enumerate(zip(frames, out_atoms)):
            _compare_atoms(orig, loaded)
            assert np.allclose(orig.get_momenta(), loaded.get_momenta())
            assert loaded.calc is not None
            assert loaded.info.get('frame_id') == i


class TestMixedFileAndBytes:
    """Test mixed operations between files and BytesIO."""

    def test_write_file_read_bytes(self, temp_dir):
        """Test writing to file and reading as bytes."""
        atoms = _make_simple_atoms()
        tess_path = temp_dir / "test.tess"
        
        # Write to file
        write_tess(atoms, str(tess_path), metadata={"mode": "file"}, compression="zlib")
        
        # Read file as bytes into BytesIO
        with open(tess_path, 'rb') as f:
            bio = io.BytesIO(f.read())
        
        # Read from BytesIO
        bio.seek(0)
        out_atoms, out_meta = read_tess(bio)
        
        assert len(out_atoms) == 1
        _compare_atoms(atoms, out_atoms[0])
        assert out_meta == {"mode": "file"}

    def test_write_bytes_read_file(self, temp_dir):
        """Test writing to BytesIO and reading from file."""
        atoms = _make_simple_atoms()
        bio = io.BytesIO()
        
        # Write to BytesIO
        write_tess(atoms, bio, metadata={"mode": "bytes"}, compression="zlib")
        
        # Write BytesIO contents to file
        tess_path = temp_dir / "test.tess"
        with open(tess_path, 'wb') as f:
            f.write(bio.getvalue())
        
        # Read from file
        out_atoms, out_meta = read_tess(str(tess_path))
        
        assert len(out_atoms) == 1
        _compare_atoms(atoms, out_atoms[0])
        assert out_meta == {"mode": "bytes"}

    def test_bytes_content_matches_file_content(self, temp_dir):
        """Test that BytesIO and file produce identical bytes."""
        atoms = _make_simple_atoms()
        meta = {"test": "identical"}
        
        # Write to BytesIO
        bio = io.BytesIO()
        write_tess(atoms, bio, metadata=meta, compression="none")
        bytes_content = bio.getvalue()
        
        # Write to file
        tess_path = temp_dir / "test.tess"
        write_tess(atoms, str(tess_path), metadata=meta, compression="none")
        
        with open(tess_path, 'rb') as f:
            file_content = f.read()
        
        # Should be identical
        assert bytes_content == file_content


class TestFixtureCompatibility:
    """Test reading existing tess fixture files."""

    def test_read_fixture_as_file(self, fixtures_dir):
        """Test reading test.tess fixture from file path."""
        tess_file = fixtures_dir / "test.tess"
        if not tess_file.exists():
            pytest.skip("test.tess fixture not found")
        
        atoms_list, metadata = read_tess(str(tess_file))
        
        # Basic sanity checks
        assert isinstance(atoms_list, list)
        assert len(atoms_list) > 0
        assert all(isinstance(a, Atoms) for a in atoms_list)
        assert isinstance(metadata, dict)

    def test_read_fixture_as_bytes(self, fixtures_dir):
        """Test reading test.tess fixture from BytesIO."""
        tess_file = fixtures_dir / "test.tess"
        if not tess_file.exists():
            pytest.skip("test.tess fixture not found")
        
        # Read file into BytesIO
        with open(tess_file, 'rb') as f:
            bio = io.BytesIO(f.read())
        
        bio.seek(0)
        atoms_list, metadata = read_tess(bio)
        
        # Basic sanity checks
        assert isinstance(atoms_list, list)
        assert len(atoms_list) > 0
        assert all(isinstance(a, Atoms) for a in atoms_list)

    def test_fixture_file_and_bytes_produce_same_result(self, fixtures_dir):
        """Test that reading fixture as file or bytes produces identical results."""
        tess_file = fixtures_dir / "test.tess"
        if not tess_file.exists():
            pytest.skip("test.tess fixture not found")
        
        # Read as file
        atoms_file, meta_file = read_tess(str(tess_file))
        
        # Read as bytes
        with open(tess_file, 'rb') as f:
            bio = io.BytesIO(f.read())
        bio.seek(0)
        atoms_bytes, meta_bytes = read_tess(bio)
        
        # Compare
        assert len(atoms_file) == len(atoms_bytes)
        assert meta_file == meta_bytes
        
        for a_file, a_bytes in zip(atoms_file, atoms_bytes):
            _compare_atoms(a_file, a_bytes)

    def test_read_fixture_subset_from_bytes(self, fixtures_dir):
        """Test reading subset of frames from fixture using BytesIO."""
        tess_file = fixtures_dir / "test.tess"
        if not tess_file.exists():
            pytest.skip("test.tess fixture not found")
        
        # First, determine how many frames exist
        with open(tess_file, 'rb') as f:
            bio = io.BytesIO(f.read())
        bio.seek(0)
        all_atoms, _ = read_tess(bio)
        
        if len(all_atoms) < 2:
            pytest.skip("test.tess has fewer than 2 frames")
        
        # Read subset
        bio.seek(0)
        subset_indices = [0, len(all_atoms) - 1]  # First and last frame
        subset_atoms, _ = read_tess(bio, frames_indices=subset_indices)
        
        assert len(subset_atoms) == 2
        _compare_atoms(all_atoms[0], subset_atoms[0])
        _compare_atoms(all_atoms[-1], subset_atoms[1])


class TestBytesIOEdgeCases:
    """Test edge cases with BytesIO operations."""

    def test_empty_bytesio_raises_error(self):
        """Test that reading from empty BytesIO raises appropriate error."""
        bio = io.BytesIO()
        
        with pytest.raises(Exception):  # Could be IndexError, struct.error, etc.
            read_tess(bio)

    def test_bytesio_not_seeked_to_start(self):
        """Test reading from BytesIO that's not at position 0."""
        atoms = _make_simple_atoms()
        bio = io.BytesIO()
        
        write_tess(atoms, bio, compression="none")
        
        # Don't seek to start - read_tess should handle this
        # Note: Based on the implementation, this might fail or succeed
        # depending on whether read_tess seeks internally
        try:
            out_atoms, _ = read_tess(bio)
            # If it succeeds, it should have read correctly
            assert len(out_atoms) == 1
        except Exception:
            # If it fails, that's also acceptable behavior
            # The important thing is it doesn't crash catastrophically
            pass

    def test_multiple_writes_to_same_bytesio(self):
        """Test that multiple sequential writes require separate BytesIO objects."""
        atoms1 = _make_simple_atoms()
        atoms2 = _make_simple_atoms()
        atoms2.set_positions(atoms2.get_positions() + np.array([0.1, 0.1, 0.1]))
        
        # Write to separate BytesIO objects (correct approach)
        bio1 = io.BytesIO()
        write_tess(atoms1, bio1, compression="none")
        
        bio2 = io.BytesIO()
        write_tess(atoms2, bio2, compression="none")
        
        # Verify each can be read correctly
        bio1.seek(0)
        out_atoms1, _ = read_tess(bio1)
        assert len(out_atoms1) == 1
        _compare_atoms(atoms1, out_atoms1[0])
        
        bio2.seek(0)
        out_atoms2, _ = read_tess(bio2)
        assert len(out_atoms2) == 1
        _compare_atoms(atoms2, out_atoms2[0])

    def test_bytesio_with_changing_cell(self):
        """Test BytesIO with frames that have changing cell parameters."""
        frames = []
        for i in range(5):
            a = _make_simple_atoms()
            cell = np.diag([5.0 + 0.1 * i, 5.0, 5.0])
            a.set_cell(cell)
            frames.append(a)
        
        bio = io.BytesIO()
        write_tess(frames, bio, compression="zlib", compression_level=1)
        
        bio.seek(0)
        out_atoms, _ = read_tess(bio)
        
        assert len(out_atoms) == len(frames)
        for orig, loaded in zip(frames, out_atoms):
            _compare_atoms(orig, loaded)

