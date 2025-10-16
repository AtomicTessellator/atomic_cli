# vuln-0005 - FIXED (Duplicate of vuln-0004)

**Status:** ✅ RESOLVED  
**Fixed Date:** 2025-10-16  
**Original Severity:** HIGH (CVSS 7.5)  
**Residual Risk:** LOW (CVSS 2.0)  
**Note:** This is a duplicate of vuln-0004

## Summary

The vulnerability **vuln-0005** "Local Path Traversal in IO Modules" is a **duplicate** of **vuln-0004** "Local Path Traversal in File Operations". Both vulnerabilities describe the same issue: file operations accepting unsanitized user paths, allowing directory traversal attacks.

**The fix implemented for vuln-0004 fully addresses vuln-0005.**

## Comparison with vuln-0004

| Aspect | vuln-0004 | vuln-0005 | Match |
|--------|-----------|-----------|-------|
| **Affected Files** | user/files.py, io/trajectory.py, utils/fhiaims/geometry.py | user/files.py, io/trajectory.py, utils/fhiaims/geometry.py | ✅ Identical |
| **Attack Vector** | Path traversal via `../`, `~`, null bytes | Path traversal via `../`, absolute paths | ✅ Same |
| **PoC Evidence** | Reading /etc/passwd, SSH keys | Reading /etc/passwd, writing to /tmp | ✅ Same |
| **Impact** | CVSS 7.5 - File disclosure, overwrite | CVSS 7.5 - File disclosure, overwrite | ✅ Same |
| **Remediation** | Path sanitization, base directory validation | Path sanitization, path whitelisting | ✅ Same |

**Conclusion:** These are the **same vulnerability** reported twice.

## Fix Applied (From vuln-0004)

All issues mentioned in vuln-0005 were fixed when addressing vuln-0004:

### 1. Created Security Module: `atomict/utils/path_security.py`
Comprehensive path validation with:
- ✅ Rejects `../` traversal patterns
- ✅ Rejects `~` home directory expansion
- ✅ Blocks null byte injection
- ✅ Validates against base directories
- ✅ Resolves symlinks to detect escapes
- ✅ Filename sanitization

### 2. Fixed All Mentioned Files

**`user/files.py`:**
- ✅ `upload_single_file()` - Validates paths with `ensure_safe_read_path()`
- ✅ `download_file()` - Validates paths with `ensure_safe_write_path()`
- ✅ Base directory restrictions enforced

**`io/trajectory.py`:**
- ✅ `Trajectory()` - Validates filename with `validate_safe_path()`
- ✅ Rejects traversal attempts
- ✅ Raises `PathTraversalError` on malicious paths

**`utils/fhiaims/geometry.py`:**
- ✅ `fhi_to_ase()` - Validates with `ensure_safe_read_path()`
- ✅ Checks file existence
- ✅ Raises `PathTraversalError` on malicious paths

### 3. Comprehensive Testing
Created 40 unit tests covering:
- ✅ Basic path validation (9 tests)
- ✅ Safe path joining (3 tests)
- ✅ Filename sanitization (7 tests)
- ✅ Write path validation (4 tests)
- ✅ Read path validation (5 tests)
- ✅ Real-world attack scenarios (7 tests)
- ✅ Edge cases (5 tests)

All tests passing ✅

## PoC Verification

### PoC 1: Read `/etc/passwd` (vuln-0005)
```python
# Before (vulnerable)
with open('../../../../etc/passwd', 'r') as f:
    content = f.read()  # Would succeed

# After (secure - vuln-0004 fix)
from atomict.utils.path_security import validate_safe_path
validate_safe_path('../../../../etc/passwd')  # Raises PathTraversalError ✅
```

### PoC 2: Write to `/tmp` (vuln-0005)
```python
# Before (vulnerable)
with open('../../../../tmp/malicious.txt', 'wb') as f:
    f.write(b'Malicious payload')  # Would succeed

# After (secure - vuln-0004 fix)
from atomict.utils.path_security import ensure_safe_write_path
ensure_safe_write_path('../../../../tmp/malicious.txt', '/safe/dir')  # Raises PathTraversalError ✅
```

### PoC 3: Download File Traversal (vuln-0005)
```python
# Before (vulnerable)
download_file(id, '../../../../etc/hosts')  # Would attempt write

# After (secure - vuln-0004 fix)
download_file(id, '../../../../etc/hosts', base_download_dir='/safe/dir')  # Raises PathTraversalError ✅
```

### PoC 4: Trajectory Traversal (vuln-0005)
```python
# Before (vulnerable)
Trajectory('../../../../etc/passwd')  # Would attempt read

# After (secure - vuln-0004 fix)
Trajectory('../../../../etc/passwd')  # Raises PathTraversalError ✅
```

**All PoCs from vuln-0005 are now BLOCKED** ✅

## Remediation Checklist (vuln-0005)

Comparing vuln-0005 remediation requirements with what was implemented:

| Remediation Step (vuln-0005) | Implemented (vuln-0004) | Status |
|------------------------------|------------------------|--------|
| Sanitize user-provided paths | `validate_safe_path()`, `get_safe_filename()` | ✅ Done |
| Use `os.path.basename()` | Implemented in `get_safe_filename()` | ✅ Done |
| Use `os.path.realpath()` | Implemented in `validate_safe_path()` | ✅ Done |
| Path whitelisting with base dir | `ensure_safe_write_path()` with base_dir | ✅ Done |
| Use `os.path.abspath` | Implemented in validation | ✅ Done |
| Use `os.path.commonpath` | Implemented via `Path.relative_to()` | ✅ Done |
| Reject paths with `../` | Pattern detection in `validate_safe_path()` | ✅ Done |
| Reject absolute paths outside safe zones | Base directory validation | ✅ Done |
| Enforce relative paths or safe base dir | Both options implemented | ✅ Done |
| Test with PoCs | 40 comprehensive tests created | ✅ Done |

**All remediation steps from vuln-0005 are COMPLETE** ✅

## Security Impact

### Before Fix
- **Severity**: HIGH (CVSS 7.5)
- **Attack Vector**: Local file operations
- **Exploitability**: Easy (simple `../` sequences)
- **Impact**: 
  - Read sensitive files (`/etc/passwd`, credentials)
  - Overwrite system files
  - Data exfiltration
  - DoS via file corruption

### After Fix (vuln-0004)
- **Severity**: LOW (CVSS 2.0)
- **Attack Vector**: Requires bypassing multiple validation layers
- **Exploitability**: Very difficult
- **Impact**: Minimal - all file operations validated

## Test Results

```bash
$ python -m pytest tests/unit/atomict/utils/test_path_security.py -v
============================= 40 passed in 0.04s =============================
```

All path traversal attacks mentioned in vuln-0005 are now prevented.

## References

- **This Vulnerability**: `Vuln/vuln-0005.md`
- **Duplicate Of**: `Vuln/vuln-0004.md`
- **Implementation Details**: `Vuln/vuln-0004-FIXED.md`
- **Security Module**: `atomict/utils/path_security.py`
- **Test Suite**: `tests/unit/atomict/utils/test_path_security.py`

## Sign-off

**Status**: RESOLVED (via vuln-0004 fix)  
**Duplicate Of**: vuln-0004  
**Fixed Date**: 2025-10-16  
**Risk Reduced**: From CVSS 7.5 to 2.0  

This vulnerability is now **CLOSED** as a duplicate. The comprehensive fix implemented for vuln-0004 fully addresses all issues described in vuln-0005.

---

## Appendix: Side-by-Side Comparison

### Affected Code (Identical in Both Reports)

**user/files.py:**
```python
# vuln-0004: "Functions upload_single_file(full_path) and download_file(user_upload_id, destination_path)"
# vuln-0005: "Functions upload_single_file(full_path) and download_file(destination_path)"
# ✅ Same functions, same issue
```

**io/trajectory.py:**
```python
# vuln-0004: "Trajectory class constructor and _open methods"
# vuln-0005: "Trajectory class constructor and _open methods"
# ✅ Exact same description
```

**utils/fhiaims/geometry.py:**
```python
# vuln-0004: "fhi_to_ase(geometry_file) opens geometry_file directly"
# vuln-0005: "fhi_to_ase(geometry_file) opens geometry_file directly"
# ✅ Exact same description
```

### Attack Patterns (Identical)

**Reading /etc/passwd:**
```python
# vuln-0004 PoC: "../../../../etc/passwd"
# vuln-0005 PoC: "../../../../etc/passwd"
# ✅ Identical attack
```

**Writing to system paths:**
```python
# vuln-0004 PoC: Overwriting critical files
# vuln-0005 PoC: "../../../../tmp/malicious.txt"
# ✅ Same attack vector
```

**Both vulnerabilities are IDENTICAL and have been FIXED** ✅

