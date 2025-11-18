import click
from pathlib import Path


@click.group()
def traj():
    """Operations on trajectory files."""
    pass


@traj.command()
@click.argument('file_path', type=click.Path(exists=True))
@click.option('-s', '--step', type=int, required=True, help='Subsample every Nth frame')
def subsample(file_path, step):
    """Subsample a trajectory file by taking every Nth frame.
    
    The output will be saved as {original_filename}_subsampled_{step}.traj
    """
    from ase.io import read, write
    from ase.io.trajectory import TrajectoryWriter

    if step <= 0:
        click.echo("Error: Step must be a positive integer")
        return

    try:
        # Parse the input file path
        input_path = Path(file_path)
        
        # Read all frames from the trajectory file
        atoms_list = read(file_path, index=":")
        
        # Subsample the frames with the given step
        subsampled_atoms = atoms_list[::step]
        
        if not subsampled_atoms:
            click.echo("Warning: No frames were selected with the given step size")
            return
        
        # Create output filename
        stem = input_path.stem
        output_file = input_path.with_name(f"{stem}_subsampled_{step}.traj")
        
        # Write the subsampled trajectory
        write(str(output_file), subsampled_atoms)
        
        click.echo(f"Successfully subsampled trajectory with step {step}")
        click.echo(f"Wrote {len(subsampled_atoms)} frames to {output_file}")
        click.echo(f"Original trajectory had {len(atoms_list)} frames")
        
    except Exception as e:
        click.echo(f"Error processing trajectory file: {str(e)}") 