"""
Secure path validation utilities to prevent path traversal attacks.

This module provides functions to safely validate and sanitize file paths
to prevent directory traversal vulnerabilities.
"""
import os
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class PathTraversalError(ValueError):
    """Raised when a path traversal attempt is detected."""
    pass


def is_safe_path(path: str, base_dir: Optional[str] = None) -> bool:
    """
    Check if a path is safe (doesn't attempt directory traversal).
    
    Args:
        path: The path to validate
        base_dir: Optional base directory to restrict paths to
        
    Returns:
        True if the path is safe, False otherwise
    """
    try:
        validate_safe_path(path, base_dir)
        return True
    except PathTraversalError:
        return False


def validate_safe_path(path: str, base_dir: Optional[str] = None) -> str:
    """
    Validate and sanitize a file path to prevent directory traversal.
    
    This function:
    1. Resolves the path to its absolute form
    2. Checks for path traversal attempts
    3. If base_dir is provided, ensures the path is within that directory
    4. Returns the validated absolute path
    
    Args:
        path: The path to validate
        base_dir: Optional base directory to restrict paths to.
                 If provided, the path must be within this directory.
        
    Returns:
        The validated absolute path
        
    Raises:
        PathTraversalError: If path traversal is detected
        ValueError: If path contains invalid characters or patterns
    """
    if not path:
        raise ValueError("Path cannot be empty")
    
    # Check for null bytes (can be used to bypass checks)
    if '\x00' in path:
        raise ValueError("Path contains null bytes")
    
    # Convert to Path object for better handling
    path_obj = Path(path)
    
    # Check for suspicious patterns
    suspicious_patterns = ['..', '~']
    path_str = str(path_obj)
    
    # Check each path component for suspicious patterns
    for part in path_obj.parts:
        if part in suspicious_patterns:
            logger.warning(f"Path traversal attempt detected: {path}")
            raise PathTraversalError(
                f"Path contains suspicious pattern '{part}': {path}"
            )
    
    # If base_dir is provided, validate against it
    if base_dir:
        # Resolve both paths to absolute form
        base_path = Path(base_dir).resolve()
        
        # If path is relative, join it with base_dir
        if not path_obj.is_absolute():
            full_path = (base_path / path_obj).resolve()
        else:
            full_path = path_obj.resolve()
        
        # Check if the resolved path is within base_dir
        try:
            full_path.relative_to(base_path)
        except ValueError:
            logger.warning(
                f"Path traversal attempt: {path} resolves outside base directory {base_dir}"
            )
            raise PathTraversalError(
                f"Path '{path}' is outside allowed directory '{base_dir}'"
            )
        
        return str(full_path)
    else:
        # No base_dir restriction, just resolve and check for traversal patterns
        resolved_path = path_obj.resolve()
        return str(resolved_path)


def safe_join(base_dir: str, *paths: str) -> str:
    """
    Safely join paths, ensuring the result stays within base_dir.
    
    This is similar to os.path.join but validates the result to prevent
    path traversal attacks.
    
    Args:
        base_dir: The base directory
        *paths: Path components to join
        
    Returns:
        The validated joined path
        
    Raises:
        PathTraversalError: If the resulting path escapes base_dir
    """
    if not paths:
        return str(Path(base_dir).resolve())
    
    # Join all paths
    joined = os.path.join(*paths)
    
    # Validate the joined path
    return validate_safe_path(joined, base_dir)


def get_safe_filename(filename: str) -> str:
    """
    Extract a safe filename from a path, removing any directory components.
    
    This extracts just the filename part and validates it doesn't contain
    suspicious characters.
    
    Args:
        filename: The filename or path
        
    Returns:
        Just the filename component, sanitized
        
    Raises:
        ValueError: If the filename is invalid
    """
    if not filename:
        raise ValueError("Filename cannot be empty")
    
    # Check for null bytes
    if '\x00' in filename:
        raise ValueError("Filename contains null bytes")
    
    # Extract just the filename (basename)
    safe_name = os.path.basename(filename)
    
    # Ensure we got something after basename extraction
    if not safe_name or safe_name in ('.', '..'):
        raise ValueError(f"Invalid filename: {filename}")
    
    # Check for suspicious characters in the filename itself
    # (These shouldn't be in a filename even without traversal)
    suspicious_chars = ['/', '\\', '\x00']
    for char in suspicious_chars:
        if char in safe_name:
            raise ValueError(f"Filename contains invalid character: {char}")
    
    return safe_name


def ensure_safe_write_path(path: str, base_dir: str, create_dirs: bool = False) -> str:
    """
    Validate a path is safe for writing and optionally create parent directories.
    
    Args:
        path: The path to validate
        base_dir: Base directory to restrict writes to
        create_dirs: If True, create parent directories if they don't exist
        
    Returns:
        The validated absolute path
        
    Raises:
        PathTraversalError: If path validation fails
        OSError: If directory creation fails
    """
    # Validate the path
    safe_path = validate_safe_path(path, base_dir)
    
    # Create parent directories if requested
    if create_dirs:
        parent_dir = os.path.dirname(safe_path)
        if parent_dir and not os.path.exists(parent_dir):
            # Validate parent directory is also within base_dir
            validate_safe_path(parent_dir, base_dir)
            os.makedirs(parent_dir, mode=0o755, exist_ok=True)
            logger.debug(f"Created directory: {parent_dir}")
    
    return safe_path


def ensure_safe_read_path(path: str, base_dir: Optional[str] = None) -> str:
    """
    Validate a path is safe for reading.
    
    Args:
        path: The path to validate
        base_dir: Optional base directory to restrict reads to
        
    Returns:
        The validated absolute path
        
    Raises:
        PathTraversalError: If path validation fails
        FileNotFoundError: If the file doesn't exist
    """
    # Validate the path
    safe_path = validate_safe_path(path, base_dir)
    
    # Check file exists
    if not os.path.exists(safe_path):
        raise FileNotFoundError(f"File not found: {path}")
    
    # Check it's actually a file (not a directory)
    if not os.path.isfile(safe_path):
        raise ValueError(f"Path is not a file: {path}")
    
    return safe_path

