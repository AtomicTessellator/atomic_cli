from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_k8s_cluster(cluster_id: str, **params) -> dict:
    return retrieve_resource("api/k8s-cluster", cluster_id, **params)


def list_k8s_clusters(**params) -> dict:
    return list_resources("api/k8s-cluster", **params)


def update_k8s_cluster(cluster_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/k8s-cluster", cluster_id, fields)


def delete_k8s_cluster(cluster_id: str) -> dict:
    return delete_resource("api/k8s-cluster", cluster_id)


def get_k8s_job(job_id: str, **params) -> dict:
    return retrieve_resource("api/k8s-job", job_id, **params)


def list_k8s_jobs(**params) -> dict:
    return list_resources("api/k8s-job", **params)


def create_k8s_job(payload: dict[str, object]) -> dict:
    return create_resource("api/k8s-job", payload)


def update_k8s_job(job_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/k8s-job", job_id, fields)


def delete_k8s_job(job_id: str) -> dict:
    return delete_resource("api/k8s-job", job_id)
