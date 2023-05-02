from atomict.datasets.public import MAPPINGS
from atomict.utils.requests import download_large_file


import os
import pandas as pd


def get_public_dataset(stub: str) -> pd.DataFrame:
    """Get a public dataset from the atomict-raw-datasets s3 bucket.
    Parameters
    ----------
    stub : str
        The name of the dataset to load.
    Returns
    -------
    pd.DataFrame
        The dataset.
    """

    if stub not in MAPPINGS.keys():
        raise ValueError(f'Unknown dataset {stub}.')

    return pd.read_csv(MAPPINGS[stub])


def download_public_dataset(stub: str, destination: str) -> None:
    """Download a public dataset from the atomict-raw-datasets s3 bucket.
    Parameters
    ----------
    stub : str
        The name of the dataset to load.
    destination : str
        The path to save the dataset to.
    Returns
    -------
    None
    """

    if stub not in MAPPINGS.keys():
        raise ValueError(f'Unknown dataset {stub}.')

    url = MAPPINGS[stub]
    full_path = os.path.join(destination, url.split('/')[-1])

    print(f'Downloading {url} to {destination}...')
    download_large_file(url, full_path)
