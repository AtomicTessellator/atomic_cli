import logging
import os
import shutil

from atomict.infra.distwork.task import TaskStatus, update_task_status
from atomict.io.utils import human_filesize
from atomict.user.files import download_file, upload_single_file


def display_name(user_upload):
    if "users_name" in user_upload and user_upload["users_name"] not in ["", None]:
        return user_upload["users_name"]
    return user_upload["orig_name"]


def clear_workspace(sim, base_path: str = "./workspace"):

    target_dir = os.path.join(base_path, sim["id"])

    if os.path.exists(target_dir):
        logging.warning(f"Removing existing workspace folder {target_dir}")
        shutil.rmtree(target_dir)


def download_workspace(workspace_files, target_directory: str):
    """Download a workspace to a target directory

    Args:
        workspace_files list of UserUpload objects: List of files to download
        target_directory (str): The directory to download the files to
    """

    os.makedirs(target_directory, exist_ok=True)

    total_bytes = sum([f["user_upload"]["size"] for f in workspace_files])
    finished_bytes = 0
    for sim_file in workspace_files:
        logging.info(
            f"Downloading file {display_name(sim_file['user_upload'])} ({human_filesize(sim_file['user_upload']['size'])})"
        )

        download_file(
            sim_file["user_upload"]["id"],
            f"{target_directory}/{sim_file['user_upload']['users_name']}",
        )

        finished_bytes += sim_file["user_upload"]["size"]
        logging.info(
            f"Downloaded {human_filesize(finished_bytes)} of {human_filesize(total_bytes)} ({finished_bytes/total_bytes*100:.1f}%)"
        )


def upload_workspace(
    sim, associate_function, workspace_folder: str, starting_percent: int = 80
):
    """
    Uploads a workspace folder to the Atomic platform, and associates the uploaded files with the given simulation.

    Args:
        sim: The simulation object to associate the uploaded files with
        associate_function: The function to associate the uploaded files with the simulation (e.g. associate_user_upload_with_qe_simulation)
        workspace_folder: The folder to upload
        starting_percent: The starting percentage to use when updating the task status
    """

    simulation_id = sim["id"]

    total_size = sum(
        os.path.getsize(os.path.join(root, file))
        for root, _, files in os.walk(workspace_folder)
        for file in files
    )
    uploaded_size = 0
    last_update_percent = starting_percent

    update_task_status(sim["task"]["id"],
        percent=last_update_percent,
        progress_indeterminate=False,
    )

    for root, _, files in os.walk(workspace_folder):
        for file in files:
            file_path = os.path.join(root, file)
            inner_workspace = file_path.replace(workspace_folder, "")
            file_size = os.path.getsize(file_path)

            try:
                result = upload_single_file(file_path, inner_workspace)
                if result["status"] != "OK":
                    logging.error(f"Failed to upload {inner_workspace}")
                    logging.error(result)
                    raise Exception(f"Failed to upload {inner_workspace} {result}")
                else:
                    logging.info(f"Uploaded {inner_workspace} OK")
            except Exception as e:
                logging.error(f"Failed to upload {inner_workspace}")
                logging.error(e)
                raise

            associate_function(result["UserUpload"]["id"], simulation_id)

            uploaded_size += file_size
            current_percent = 80 + int((uploaded_size / total_size) * 20)

            # Update status if at least 2% has changed
            if current_percent - last_update_percent >= 2:
                update_task_status(sim["task"]["id"], percent=current_percent)
                last_update_percent = current_percent
