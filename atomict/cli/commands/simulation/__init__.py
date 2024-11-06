# cli/commands/simulation/__init__.py
import click

from . import catalysis, fhiaims  # qe,; phonon,; kpoint,


@click.group(name="simulation")
def simulation_group():
    """Manage all types of simulations"""
    pass


# Add simulation type subcommands
simulation_group.add_command(catalysis.catalysis)
# simulation_group.add_command(qe.qe)
simulation_group.add_command(fhiaims.fhiaims_group)
# simulation_group.add_command(phonon.phonon)
# simulation_group.add_command(kpoint.kpoint)

simulation = simulation_group
