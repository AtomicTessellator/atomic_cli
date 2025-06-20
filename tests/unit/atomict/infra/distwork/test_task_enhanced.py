import pytest
from unittest.mock import patch, Mock
from atomict.infra.distwork.task import (
    cancel_task,
    get_task_status_history,
    tail_task_logs,
    TaskStatus
)


class TestEnhancedTaskManagement:
    """Test enhanced task management functions."""
    
    @patch('atomict.infra.distwork.task.patch')
    def test_cancel_task(self, mock_patch):
        """Test task cancellation."""
        mock_patch.return_value = {
            "id": "task-123",
            "user_aborted_flag": True,
            "status": TaskStatus.USER_ABORTED.value
        }
        
        result = cancel_task("task-123")
        
        expected_payload = {"user_aborted_flag": True}
        mock_patch.assert_called_once_with("api/tasks/task-123/", payload=expected_payload)
        assert result["user_aborted_flag"] is True

    @patch('atomict.infra.distwork.task.get')
    def test_get_task_status_history_with_results(self, mock_get):
        """Test task status history retrieval with results."""
        mock_get.return_value = {
            "results": [
                {
                    "id": "history-1",
                    "task": "task-123",
                    "new_status": TaskStatus.READY.value,
                    "timestamp": "2024-01-01T00:00:00Z"
                },
                {
                    "id": "history-2", 
                    "task": "task-123",
                    "new_status": TaskStatus.RUNNING.value,
                    "timestamp": "2024-01-01T01:00:00Z"
                }
            ]
        }
        
        result = get_task_status_history("task-123")
        
        mock_get.assert_called_once_with("api/task-status-history/?task__id=task-123")
        assert len(result) == 2
        assert result[0]["new_status"] == TaskStatus.READY.value
        assert result[1]["new_status"] == TaskStatus.RUNNING.value

    @patch('atomict.infra.distwork.task.get')
    def test_get_task_status_history_direct_response(self, mock_get):
        """Test task status history when API returns direct response."""
        mock_get.return_value = [
            {"id": "history-1", "new_status": TaskStatus.COMPLETED.value}
        ]
        
        result = get_task_status_history("task-123")
        
        mock_get.assert_called_once_with("api/task-status-history/?task__id=task-123")
        assert isinstance(result, list)
        assert result[0]["new_status"] == TaskStatus.COMPLETED.value

    @patch('atomict.infra.distwork.task.get')
    def test_get_task_status_history_empty_results(self, mock_get):
        """Test task status history with empty results."""
        mock_get.return_value = {"results": []}
        
        result = get_task_status_history("task-123")
        
        assert result == []

    @patch('atomict.infra.distwork.task.get')
    def test_tail_task_logs(self, mock_get):
        """Test task log retrieval."""
        mock_get.return_value = {
            "logs": "2024-01-01 00:00:00 - Task started\n2024-01-01 00:01:00 - Processing data\n2024-01-01 00:02:00 - Task completed"
        }
        
        result = tail_task_logs("task-123")
        
        mock_get.assert_called_once_with("api/task/task-123/logs/")
        assert "logs" in result
        assert "Task started" in result["logs"]

    @patch('atomict.infra.distwork.task.get')
    def test_tail_task_logs_empty(self, mock_get):
        """Test task log retrieval with empty logs."""
        mock_get.return_value = {"logs": ""}
        
        result = tail_task_logs("task-123")
        
        assert result["logs"] == ""

    @patch('atomict.infra.distwork.task.get')
    def test_tail_task_logs_no_logs_field(self, mock_get):
        """Test task log retrieval when logs field is missing."""
        mock_get.return_value = {"message": "No logs available"}
        
        result = tail_task_logs("task-123")
        
        assert "message" in result
        assert result["message"] == "No logs available"
