from enum import Enum

from atomict.api import get, patch
from atomict.resource_helpers import list_resources, update_resource
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


def get_task(task_id: str):
    return get(f"api/tasks/{task_id}/")


def list_tasks(**params):
    if "depth" not in params:
        params["depth"] = 2
    return list_resources("api/tasks", **params)


def get_task_status_history(task_id: str, **params):
    query_params = {"id": task_id}
    query_params.update(params)
    return list_resources("api/task-status-history", **query_params)


def cancel_task(task_id: str):
    return update_resource("api/tasks", task_id, {"status": TaskStatus.USER_ABORTED.value})


def task_should_abort(task_id: str) -> bool:
    task = get_task(task_id)
    return task["user_aborted_flag"]


def except_on_user_abort(task_id: str):
    if task_should_abort(task_id):
        raise UserTaskAbortException(f"User aborted task {task_id}")


def update_task_status(
    task_id: str, status: TaskStatus = None, error_msg: str = None, percent: int = None, progress_indeterminate: bool = None
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

    res = patch(f"api/tasks/{task_id}/", payload=payload)
    return res
