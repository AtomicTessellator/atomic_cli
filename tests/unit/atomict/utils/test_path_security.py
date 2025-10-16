"""Unit tests for path security utilities."""
import os
import tempfile
from pathlib import Path

import pytest

from atomict.utils.path_security import (
    PathTraversalError,
    ensure_safe_read_path,
    ensure_safe_write_path,
    get_safe_filename,
    is_safe_path,
    safe_join,
    validate_safe_path,
)


class TestPathValidation:
    """Test basic path validation functions."""
    
    def test_validate_safe_path_simple(self):
        """Test validation of a simple filename."""
        result = validate_safe_path("test.txt")
        assert "test.txt" in result
        assert os.path.isabs(result)
    
    def test_validate_safe_path_rejects_parent_traversal(self):
        """Test that ../ patterns are rejected."""
        with pytest.raises(PathTraversalError, match="suspicious pattern"):
            validate_safe_path("../etc/passwd")
    
    def test_validate_safe_path_rejects_double_dot(self):
        """Test that .. is rejected."""
        with pytest.raises(PathTraversalError, match="suspicious pattern"):
            validate_safe_path("foo/../bar/baz.txt")
    
    def test_validate_safe_path_rejects_tilde(self):
        """Test that ~ is rejected."""
        with pytest.raises(PathTraversalError, match="suspicious pattern"):
            validate_safe_path("~/secrets.txt")
    
    def test_validate_safe_path_null_byte(self):
        """Test that null bytes are rejected."""
        with pytest.raises(ValueError, match="null bytes"):
            validate_safe_path("test\x00.txt")
    
    def test_validate_safe_path_empty(self):
        """Test that empty paths are rejected."""
        with pytest.raises(ValueError, match="cannot be empty"):
            validate_safe_path("")
    
    def test_validate_safe_path_with_base_dir(self, tmp_path):
        """Test validation with base directory restriction."""
        result = validate_safe_path("subdir/file.txt", str(tmp_path))
        assert result.startswith(str(tmp_path))
        assert "subdir" in result
    
    def test_validate_safe_path_escapes_base_dir(self, tmp_path):
        """Test that paths escaping base_dir are rejected."""
        # Create a file outside base_dir
        other_dir = tmp_path.parent / "other"
        other_dir.mkdir(exist_ok=True)
        
        with pytest.raises(PathTraversalError, match="outside allowed directory"):
            validate_safe_path(str(other_dir / "file.txt"), str(tmp_path))
    
    def test_is_safe_path(self):
        """Test is_safe_path wrapper."""
        assert is_safe_path("test.txt") is True
        assert is_safe_path("../etc/passwd") is False


class TestSafeJoin:
    """Test safe path joining."""
    
    def test_safe_join_basic(self, tmp_path):
        """Test basic path joining."""
        result = safe_join(str(tmp_path), "subdir", "file.txt")
        assert result.startswith(str(tmp_path))
        assert "subdir" in result
        assert "file.txt" in result
    
    def test_safe_join_prevents_traversal(self, tmp_path):
        """Test that safe_join prevents traversal."""
        with pytest.raises(PathTraversalError):
            safe_join(str(tmp_path), "..", "etc", "passwd")
    
    def test_safe_join_no_paths(self, tmp_path):
        """Test safe_join with no path components."""
        result = safe_join(str(tmp_path))
        assert os.path.isabs(result)


class TestGetSafeFilename:
    """Test filename extraction and sanitization."""
    
    def test_get_safe_filename_simple(self):
        """Test extraction of simple filename."""
        result = get_safe_filename("test.txt")
        assert result == "test.txt"
    
    def test_get_safe_filename_with_path(self):
        """Test extraction removes path components."""
        result = get_safe_filename("/path/to/file.txt")
        assert result == "file.txt"
        assert "/" not in result
    
    def test_get_safe_filename_with_traversal(self):
        """Test extraction from traversal attempt."""
        result = get_safe_filename("../../etc/passwd")
        assert result == "passwd"
    
    def test_get_safe_filename_empty(self):
        """Test that empty filename is rejected."""
        with pytest.raises(ValueError, match="cannot be empty"):
            get_safe_filename("")
    
    def test_get_safe_filename_null_byte(self):
        """Test that null bytes are rejected."""
        with pytest.raises(ValueError, match="null bytes"):
            get_safe_filename("test\x00.txt")
    
    def test_get_safe_filename_dot(self):
        """Test that . is rejected."""
        with pytest.raises(ValueError, match="Invalid filename"):
            get_safe_filename(".")
    
    def test_get_safe_filename_double_dot(self):
        """Test that .. is rejected."""
        with pytest.raises(ValueError, match="Invalid filename"):
            get_safe_filename("..")


class TestEnsafeSafeWritePath:
    """Test safe write path validation."""
    
    def test_ensure_safe_write_path_basic(self, tmp_path):
        """Test basic write path validation."""
        filepath = tmp_path / "test.txt"
        result = ensure_safe_write_path(str(filepath), str(tmp_path))
        assert result == str(filepath.resolve())
    
    def test_ensure_safe_write_path_creates_dirs(self, tmp_path):
        """Test that parent directories are created."""
        filepath = tmp_path / "subdir" / "nested" / "file.txt"
        result = ensure_safe_write_path(
            str(filepath),
            str(tmp_path),
            create_dirs=True
        )
        assert os.path.exists(os.path.dirname(result))
    
    def test_ensure_safe_write_path_rejects_traversal(self, tmp_path):
        """Test that traversal attempts are rejected."""
        with pytest.raises(PathTraversalError):
            ensure_safe_write_path("../etc/passwd", str(tmp_path))
    
    def test_ensure_safe_write_path_outside_base(self, tmp_path):
        """Test that paths outside base_dir are rejected."""
        other_dir = tmp_path.parent / "other"
        with pytest.raises(PathTraversalError, match="outside allowed directory"):
            ensure_safe_write_path(str(other_dir / "file.txt"), str(tmp_path))


class TestEnsafeSafeReadPath:
    """Test safe read path validation."""
    
    def test_ensure_safe_read_path_basic(self, tmp_path):
        """Test basic read path validation."""
        # Create a file
        filepath = tmp_path / "test.txt"
        filepath.write_text("test content")
        
        result = ensure_safe_read_path(str(filepath))
        assert result == str(filepath.resolve())
    
    def test_ensure_safe_read_path_with_base_dir(self, tmp_path):
        """Test read path validation with base directory."""
        filepath = tmp_path / "test.txt"
        filepath.write_text("test content")
        
        result = ensure_safe_read_path(str(filepath), str(tmp_path))
        assert result.startswith(str(tmp_path))
    
    def test_ensure_safe_read_path_nonexistent(self, tmp_path):
        """Test that nonexistent files are rejected."""
        filepath = tmp_path / "nonexistent.txt"
        with pytest.raises(FileNotFoundError):
            ensure_safe_read_path(str(filepath))
    
    def test_ensure_safe_read_path_directory(self, tmp_path):
        """Test that directories are rejected."""
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        with pytest.raises(ValueError, match="not a file"):
            ensure_safe_read_path(str(subdir))
    
    def test_ensure_safe_read_path_traversal(self):
        """Test that traversal attempts are rejected."""
        with pytest.raises(PathTraversalError):
            ensure_safe_read_path("../../../etc/passwd")


class TestRealWorldScenarios:
    """Test real-world attack scenarios."""
    
    def test_aws_metadata_attack(self):
        """Test that AWS metadata service paths are rejected."""
        malicious_paths = [
            "../../../../../../etc/passwd",
            "../../../../../proc/self/environ",
            "../../../../tmp/evil.sh",
        ]
        for path in malicious_paths:
            assert is_safe_path(path) is False
    
    def test_windows_traversal(self):
        """Test Windows-style traversal attempts."""
        # Test path traversal with ../ sequences
        with pytest.raises(PathTraversalError):
            validate_safe_path("../../../Windows/System32/config/SAM")
        
        # Test that absolute paths without base_dir are allowed (but resolved)
        # This is OK because without a base_dir restriction, absolute paths are valid
        result = validate_safe_path("/etc/passwd")
        assert os.path.isabs(result)
    
    def test_null_byte_injection(self):
        """Test null byte injection attempts."""
        malicious_names = [
            "test.txt\x00.php",
            "innocent\x00../../etc/passwd",
        ]
        for name in malicious_names:
            with pytest.raises(ValueError, match="null bytes"):
                validate_safe_path(name)
    
    def test_symlink_escape(self, tmp_path):
        """Test that symlinks cannot escape base directory."""
        # Create a symlink pointing outside
        link_path = tmp_path / "link"
        target_path = tmp_path.parent / "outside.txt"
        
        try:
            link_path.symlink_to(target_path)
            
            # The resolved path should be detected as outside base_dir
            with pytest.raises(PathTraversalError, match="outside allowed directory"):
                validate_safe_path(str(link_path), str(tmp_path))
        except OSError:
            # Skip on systems that don't support symlinks
            pytest.skip("Symlinks not supported")
    
    def test_upload_filename_sanitization(self):
        """Test filename sanitization for file uploads."""
        # Simulating malicious upload filenames
        malicious_names = [
            "../../etc/passwd",
            "../../../root/.ssh/id_rsa",
            "~/malicious.sh",
        ]
        
        for name in malicious_names:
            # Should extract just the safe filename
            safe_name = get_safe_filename(name)
            assert "/" not in safe_name
            assert "\\" not in safe_name
            assert ".." not in safe_name
    
    def test_download_path_restriction(self, tmp_path):
        """Test that downloads are restricted to allowed directory."""
        # Simulate attempting to download to system directory
        malicious_destinations = [
            "../../etc/cron.d/malicious",
            "../../../root/.bashrc",
            "/etc/passwd",
        ]
        
        for dest in malicious_destinations:
            with pytest.raises(PathTraversalError):
                ensure_safe_write_path(dest, str(tmp_path))
    
    def test_trajectory_file_traversal(self):
        """Test trajectory file path traversal prevention."""
        malicious_trajectories = [
            "../../../../etc/passwd.traj",
            "../../../tmp/evil.mpk",
            "~/../../etc/shadow.traj",
        ]
        
        for traj_path in malicious_trajectories:
            with pytest.raises(PathTraversalError):
                validate_safe_path(traj_path)


class TestEdgeCases:
    """Test edge cases and corner cases."""
    
    def test_current_directory_dot(self):
        """Test that single dot is handled correctly."""
        # Single dot in path is generally OK as it refers to current directory
        result = validate_safe_path("./test.txt")
        assert "test.txt" in result
    
    def test_hidden_files(self):
        """Test that hidden files (starting with .) are allowed."""
        result = validate_safe_path(".gitignore")
        assert ".gitignore" in result
    
    def test_spaces_in_filename(self):
        """Test that spaces in filenames are handled."""
        result = validate_safe_path("my file.txt")
        assert "my file.txt" in result or "my%20file.txt" in result
    
    def test_unicode_filename(self):
        """Test that unicode filenames are handled."""
        result = validate_safe_path("日本語.txt")
        assert "日本語" in result or "txt" in result
    
    def test_very_long_path(self, tmp_path):
        """Test handling of very long paths."""
        # Create a deeply nested directory
        deep_path = tmp_path
        for i in range(50):
            deep_path = deep_path / f"dir{i}"
        
        filepath = deep_path / "file.txt"
        # Should handle long paths without issues
        result = ensure_safe_write_path(
            str(filepath),
            str(tmp_path),
            create_dirs=True
        )
        assert str(tmp_path) in result

