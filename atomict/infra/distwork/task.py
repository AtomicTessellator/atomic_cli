from enum import Enum
from typing import Any, Dict, Iterator, List, Optional

from atomict.api import get, patch
from atomict.exceptions import UserTaskAbortException


class TaskStatus(Enum):
    DRAFT = 0
    READY = 1
    RUNNING = 2
    COMPLETED = 3
    ERROR = 4
    PAUSED = 5
    USER_ABORTED = 6


def get_task(task_uuid: str) -> Any:
    """Get details of a specific task."""
    return get(f"api/tasks/{task_uuid}/")


def task_should_abort(task_uuid: str) -> bool:
    """Check if a task should abort based on user flag."""
    task = get_task(task_uuid)
    return task.get("user_aborted_flag", False)


def except_on_user_abort(task_uuid: str) -> None:
    """Raise exception if task should abort."""
    if task_should_abort(task_uuid):
        raise UserTaskAbortException(f"User aborted task {task_uuid}")


def update_task_status(
    task_uuid: str,
    status: Optional[TaskStatus] = None,
    error_msg: Optional[str] = None,
    percent: Optional[int] = None,
    progress_indeterminate: Optional[bool] = None,
) -> Any:
    """Update task status and progress information."""
    payload = {}

    if status:
        payload["status"] = status.value

    if error_msg:
        payload["error"] = error_msg

    if percent:
        payload["progress"] = percent

    if progress_indeterminate is not None:
        payload["progress_indeterminate"] = progress_indeterminate

    res = patch(f"api/tasks/{task_uuid}/", payload=payload)
    return res


def cancel_task(task_uuid: str) -> Any:
    """
    Cancel a running task by setting user_aborted_flag.

    Args:
        task_uuid: UUID of the task to cancel

    Returns:
        Updated task dictionary
    """
    payload = {"user_aborted_flag": True}
    return patch(f"api/tasks/{task_uuid}/", payload=payload)


def get_task_status_history(task_uuid: str) -> Any:
    """
    Get status history for a specific task.

    Args:
        task_uuid: UUID of the task

    Returns:
        List of status history entries
    """
    response = get(f"api/task-status-history/?task__id={task_uuid}")
    if isinstance(response, dict) and "results" in response:
        return response["results"]
    return response


def tail_task_logs(task_uuid: str) -> Any:
    """
    Get logs for a specific task.

    Note: This is an alias for get_task_logs from the k8s module for convenience.
    For streaming logs, use the CLI or implement streaming in the client.

    Args:
        task_uuid: UUID of the task

    Returns:
        Dictionary containing task logs
    """
    return get(f"api/task/{task_uuid}/logs/")
