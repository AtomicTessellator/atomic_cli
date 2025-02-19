# cli/commands/ngatlas.py
from typing import Optional

import click
from rich.console import Console

from atomict.cli.commands.common import create_table
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info


console = Console()

@click.group(name="ngatlas")
def ngatlas():
    """Get NGAtlas isotope data."""
    pass


@ngatlas.command()
@click.option("--atomic-number", type=int, help="Filter by atomic number")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get(atomic_number: Optional[int] = None, json_output: bool = False, fetch_all: bool = False):
    """Get isotopes from the ngatlas"""
    params = {}
    if atomic_number is not None:
        params["atomic_number"] = atomic_number

    client = get_client()
    
    if fetch_all:
        results = client.get_all("/api/ngatlas/", params)
    else:
        results = client.get("/api/ngatlas/", params)

    if json_output:
        console.print_json(data=results)
        return

    columns = [
        ("ID", "id", None),
        ("Name", "element", lambda x: x.get("name") if x else "N/A"),
        # ("Element Category", "element", lambda x: x.get('category') if x else 'N/A'),
        # ("Element Density", "element", lambda x: str(x.get('density')) if x and x.get('density') is not None else 'N/A'),
        # ("Element Melt", "element", lambda x: str(x.get('melt')) if x and x.get('melt') is not None else 'N/A'),
        # ("Element Phase", "element", lambda x: x.get('phase') if x else 'N/A'),
        # ("Element Electron Configuration", "element", lambda x: x.get('electron_configuration') if x else 'N/A'),
        ("Filename", "fname", None),
        ("Reaction Code", "reacode", None),
    ]
    items, footer_string = get_pagination_info(results)

    table = create_table(
        columns=columns, items=items, title="Isotopes in ngatlas", caption=footer_string
    )

    console.print(table)
