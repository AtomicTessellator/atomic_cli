from enum import Enum

from atomict.api import get, patch
from atomict.exceptions import UserTaskAbortException


class SimulationAction(Enum):
    SAVE_DRAFT = "DRAFT"
    LAUNCH = "LAUNCH"


class TaskStatus(Enum):
    DRAFT = 0
    READY = 1
    RUNNING = 2
    COMPLETED = 3
    ERROR = 4
    PAUSED = 5
    USER_ABORTED = 6


def get_task(task_id: str, *, api_root: str = None, token: str = None):
    return get(f"api/tasks/{task_id}/", api_root=api_root, token=token)


def task_should_abort(task_id: str, *, api_root: str = None, token: str = None) -> bool:
    task = get_task(task_id, api_root=api_root, token=token)
    return task["user_aborted_flag"]


def except_on_user_abort(task_id: str, *, api_root: str = None, token: str = None):
    if task_should_abort(task_id, api_root=api_root, token=token):
        raise UserTaskAbortException(f"User aborted task {task_id}")


def update_task_status(
    task_id: str,
    status: TaskStatus = None,
    error_msg: str = None,
    percent: int = None,
    progress_indeterminate: bool = None,
    *,
    api_root: str = None,
    token: str = None,
):
    payload = {}

    if status:
        payload["status"] = status.value

    if error_msg:
        payload["error"] = error_msg

    if percent:
        payload["progress"] = percent

    if progress_indeterminate is not None:
        payload["progress_indeterminate"] = progress_indeterminate

    res = patch(f"api/tasks/{task_id}/", payload=payload, api_root=api_root, token=token)
    return res
