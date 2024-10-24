import click
from . import ea


@click.group(name='exploration')
def exploration_group():
    """Manage exploration tasks"""
    pass


exploration_group.add_command(ea.ea_group)
