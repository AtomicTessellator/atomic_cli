import pytest
from unittest.mock import Mock, patch, MagicMock
from atomict.user.files import (
    delete_user_upload,
    get_user_uploads,
    get_user_upload,
    update_user_upload,
    upload_multiple_files,
    get_file_upload_filter
)
from atomict.user.workspace import (
    delete_user_workspace,
    get_workspace_summary,
    clean_workspace
)


class TestUserFileManagement:
    """Test user file management functions"""
    
    @patch('atomict.user.files.delete')
    def test_delete_user_upload(self, mock_delete):
        """Test delete_user_upload function"""
        # Mock successful deletion
        mock_delete.return_value = {"success": True}
        
        result = delete_user_upload("test-uuid")
        assert result == True
        mock_delete.assert_called_once_with("user/file_upload_delete/test-uuid/")
        
        # Test when delete returns None
        mock_delete.return_value = None
        result = delete_user_upload("test-uuid")
        assert result == False
    
    @patch('atomict.user.files.get')
    def test_get_user_uploads(self, mock_get):
        """Test get_user_uploads function"""
        # Mock API response
        mock_response = {
            "results": [
                {"id": "1", "orig_name": "file1.txt"},
                {"id": "2", "orig_name": "file2.txt"}
            ],
            "count": 2
        }
        mock_get.return_value = mock_response
        
        result = get_user_uploads(search="test", limit=10)
        assert result == mock_response
        mock_get.assert_called_once_with("api/user-upload/?search=test&limit=10")
        
        # Test with no parameters
        mock_get.reset_mock()
        get_user_uploads()
        mock_get.assert_called_once_with("api/user-upload/")
    
    @patch('atomict.user.files.get')
    def test_get_user_upload(self, mock_get):
        """Test get_user_upload function"""
        mock_response = {"id": "test-uuid", "orig_name": "test.txt"}
        mock_get.return_value = mock_response
        
        result = get_user_upload("test-uuid", include_content=True)
        assert result == mock_response
        mock_get.assert_called_once_with("api/user-upload/test-uuid/?include_content=true")
    
    @patch('atomict.user.files.patch')
    def test_update_user_upload(self, mock_patch):
        """Test update_user_upload function"""
        mock_response = {"id": "test-uuid", "users_description": "Updated"}
        mock_patch.return_value = mock_response
        
        result = update_user_upload("test-uuid", users_description="Updated")
        assert result == mock_response
        mock_patch.assert_called_once_with("api/user-upload/test-uuid/", payload={"users_description": "Updated"})
    
    @patch('atomict.user.files.upload_single_file')
    def test_upload_multiple_files(self, mock_upload):
        """Test upload_multiple_files function"""
        # Mock successful uploads
        mock_upload.return_value = {"status": "OK", "UserUpload": {"id": "123"}}
        
        result = upload_multiple_files(["/path/to/file1.txt", "/path/to/file2.txt"])
        assert len(result) == 2
        assert all(item["success"] for item in result)
        assert mock_upload.call_count == 2
        
        # Test with upload failure
        mock_upload.side_effect = Exception("Upload failed")
        result = upload_multiple_files(["/path/to/file1.txt"])
        assert len(result) == 1
        assert result[0]["success"] == False
        assert "Upload failed" in result[0]["error"]


class TestWorkspaceManagement:
    """Test workspace management functions"""
    
    @patch('atomict.user.workspace.get_user_uploads')
    @patch('atomict.user.workspace.delete_user_upload')
    def test_delete_user_workspace(self, mock_delete, mock_get_uploads):
        """Test delete_user_workspace function"""
        # Mock uploads response
        mock_uploads = [
            {"uuid": "uuid1", "orig_name": "file1.txt"},
            {"id": "id2", "orig_name": "file2.txt"}  # Test fallback to id
        ]
        mock_get_uploads.return_value = {"results": mock_uploads}
        mock_delete.return_value = True
        
        result = delete_user_workspace()
        
        assert result["deleted_count"] == 2
        assert result["failed_count"] == 0
        assert len(result["errors"]) == 0
        assert mock_delete.call_count == 2
    
    @patch('atomict.user.workspace.get_user_uploads')
    def test_get_workspace_summary(self, mock_get_uploads):
        """Test get_workspace_summary function"""
        # Mock uploads response
        mock_uploads = [
            {"size": 1000, "type": "text", "orig_name": "file1.txt"},
            {"size": 2000, "type": "binary", "orig_name": "file2.bin"},
            {"size": 500, "type": "text", "orig_name": "file3.txt"}
        ]
        mock_get_uploads.return_value = {"results": mock_uploads}
        
        result = get_workspace_summary()
        
        assert result["total_files"] == 3
        assert result["total_size"] == 3500
        assert result["file_types"]["text"] == 2
        assert result["file_types"]["binary"] == 1
        assert len(result["recent_files"]) == 3
    
    @patch('atomict.user.workspace.get_user_uploads')
    @patch('atomict.user.workspace.delete_user_upload')
    def test_clean_workspace(self, mock_delete, mock_get_uploads):
        """Test clean_workspace function"""
        from datetime import datetime, timedelta
        
        # Create old and new files
        old_date = (datetime.now() - timedelta(days=40)).isoformat()
        new_date = (datetime.now() - timedelta(days=10)).isoformat()
        
        mock_uploads = [
            {"uuid": "old1", "uploaded": old_date, "type": "text", "orig_name": "old1.txt"},
            {"uuid": "new1", "uploaded": new_date, "type": "text", "orig_name": "new1.txt"},
            {"uuid": "old2", "uploaded": old_date, "type": "binary", "orig_name": "old2.bin"}
        ]
        mock_get_uploads.return_value = {"results": mock_uploads}
        mock_delete.return_value = True
        
        # Clean files older than 30 days, only text files
        result = clean_workspace(older_than_days=30, file_types=["text"])
        
        assert result["deleted_count"] == 1  # Only old1.txt should be deleted
        assert result["skipped_count"] == 2  # new1.txt (too new) and old2.bin (wrong type)
        assert result["failed_count"] == 0


class TestErrorHandling:
    """Test error handling in user management functions"""
    
    @patch('atomict.user.files.get')
    def test_get_user_uploads_none_response(self, mock_get):
        """Test handling of None response from API"""
        mock_get.return_value = None
        
        result = get_user_uploads()
        assert result == {}
    
    @patch('atomict.user.workspace.get_user_uploads')
    def test_delete_user_workspace_api_error(self, mock_get_uploads):
        """Test handling of API errors in workspace deletion"""
        mock_get_uploads.side_effect = Exception("API Error")
        
        result = delete_user_workspace()
        
        assert result["deleted_count"] == 0
        assert result["failed_count"] == 0
        assert len(result["errors"]) == 1
        assert "Failed to retrieve workspace files" in result["errors"][0]
