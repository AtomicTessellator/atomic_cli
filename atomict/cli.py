from atomict.datasets import download_public_dataset
from atomict.datasets.public import MAPPINGS


import fire
import os
from prettytable import PrettyTable
from typing import Optional


class AtomicT:
  """Atomic Tessellator CLI."""

  def list_datasets(self):
    """List all public datasets available for download."""
    table = PrettyTable(['Dataset Ident', 'Fetch From'])
    table.align = 'l'

    for ident, url in MAPPINGS.items():
      table.add_row([ident, url])

    print(table)
    print()
    print('Use atomict get <dataset_ident> to fetch a dataset.')
    print()

  def download_dataset(self, identifier: str, destination: Optional[str] = None):
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

    # If destination is supplied it must be a valid folder path
    if destination:
      if not os.path.isdir(destination):
        raise ValueError(f'Destination: {destination} must be a valid folder')
    else:
      destination = os.getcwd()

    try:
      download_public_dataset(identifier, destination)
    except ValueError as e:
      return e


def main():
    fire.Fire(AtomicT, name='atomict')


if __name__ == '__main__':
    main()

