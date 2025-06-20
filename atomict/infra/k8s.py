"""
Kubernetes cluster management functionality.

This module provides functions for viewing Kubernetes clusters
through the Atomic Tessellator API.
"""

from typing import Optional, Any
from atomict.api import get


def list_clusters(
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    **filters: Any
) -> Any:
    """
    List Kubernetes clusters.
    
    Args:
        search: Search term to filter clusters by name
        ordering: Field to order results by (e.g., 'created_at', '-created_at')
        **filters: Additional filter parameters (e.g., public=True)
        
    Returns:
        List of cluster dictionaries
    """
    params = []
    
    if search:
        params.append(f"search={search}")
    if ordering:
        params.append(f"ordering={ordering}")
    
    # Add any additional filters
    for key, value in filters.items():
        params.append(f"{key}={value}")
    
    query_string = "&".join(params)
    path = f"api/k8s-cluster/{'?' + query_string if query_string else ''}"
    
    response = get(path)
    if isinstance(response, dict) and "results" in response:
        return response["results"]
    return response


def get_cluster(cluster_id: str) -> Any:
    """
    Get details of a specific Kubernetes cluster.
    
    Args:
        cluster_id: UUID of the cluster
        
    Returns:
        Cluster details dictionary
    """
    return get(f"api/k8s-cluster/{cluster_id}/")
