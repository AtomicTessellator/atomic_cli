from atomict.datasets.io import get_public_dataset, download_public_dataset
from atomict.datasets.public import MAPPINGS

import pandas as pd
from typing import List

# Helper function so that the API looks cleaner:
#
#  atomict.datasets.get('bace')
#
# instead of:
#
# atomict.datasets.io.get_public_dataset('bace')
def get(slug: str) -> pd.DataFrame:
    """Get a public dataset from the atomict-raw-datasets s3 bucket.
    Parameters
    ----------
    slug : str
        The name of the dataset to load.
    Returns
    -------
    pd.DataFrame
        The dataset.
    """
    return get_public_dataset(slug)
