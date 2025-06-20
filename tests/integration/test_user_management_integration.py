"""
Integration tests for user management functionality.
Requires valid authentication and test environment.
"""
import os
import pytest
from atomict.user.files import (
    get_user_uploads,
    get_user_upload,
    delete_user_upload,
    upload_multiple_files
)
from atomict.user.workspace import (
    get_workspace_summary,
    delete_user_workspace,
    clean_workspace
)


@pytest.mark.integration
class TestUserManagementIntegration:
    """Integration tests for user management functions"""
    
    def test_get_user_uploads(self):
        """Test listing user uploads"""
        result = get_user_uploads(limit=5)
        
        # Result should be a dict or list
        assert result is not None
        
        # If it's a paginated response, check structure
        if isinstance(result, dict):
            assert "results" in result or "count" in result or len(result) >= 0
        elif isinstance(result, list):
            assert len(result) >= 0
    
    def test_get_workspace_summary(self):
        """Test workspace summary functionality"""
        summary = get_workspace_summary()
        
        assert isinstance(summary, dict)
        assert "total_files" in summary
        assert "total_size" in summary
        assert "file_types" in summary
        assert "recent_files" in summary
        
        # Validate data types
        assert isinstance(summary["total_files"], int)
        assert isinstance(summary["total_size"], int)
        assert isinstance(summary["file_types"], dict)
        assert isinstance(summary["recent_files"], list)
    
    def test_clean_workspace_dry_run(self):
        """Test workspace cleaning with safe parameters"""
        # Use very restrictive parameters to avoid deleting real data
        result = clean_workspace(
            older_than_days=365,  # Only very old files
            max_files=1,          # Maximum 1 file
            file_types=["nonexistent_type"]  # Non-existent type
        )
        
        assert isinstance(result, dict)
        assert "deleted_count" in result
        assert "failed_count" in result  
        assert "skipped_count" in result
        assert "errors" in result
        
        # Should skip all files due to restrictive filter
        assert result["deleted_count"] == 0
        assert isinstance(result["errors"], list)
    
    @pytest.mark.skip(reason="Requires specific upload for testing")
    def test_get_user_upload_details(self):
        """Test getting specific upload details - requires existing upload"""
        # This would need a known upload UUID to test
        # result = get_user_upload("test-uuid")
        # assert result is not None 
        pass
    
    @pytest.mark.skip(reason="Destructive operation - only run with test data")
    def test_delete_user_workspace_safe(self):
        """Test workspace deletion - DESTRUCTIVE, only run with test data"""
        # This is a destructive operation that should only be run
        # in a test environment with disposable data
        # result = delete_user_workspace()
        # assert "deleted_count" in result
        pass


if __name__ == "__main__":
    # Quick test to verify basic functionality
    print("Testing user uploads listing...")
    try:
        uploads = get_user_uploads(limit=1)
        print(f"✓ get_user_uploads() returned: {type(uploads)}")
    except Exception as e:
        print(f"✗ get_user_uploads() failed: {e}")
    
    print("\nTesting workspace summary...")
    try:
        summary = get_workspace_summary()
        print(f"✓ get_workspace_summary() returned: {summary}")
    except Exception as e:
        print(f"✗ get_workspace_summary() failed: {e}")
    
    print("\nTesting safe workspace clean...")
    try:
        result = clean_workspace(older_than_days=1000, max_files=0)
        print(f"✓ clean_workspace() returned: {result}")
    except Exception as e:
        print(f"✗ clean_workspace() failed: {e}")
