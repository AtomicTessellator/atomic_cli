#! /bin/bash

# Check if filename argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename>"
    echo "Example: $0 litacu_100 (will look for litacu_100.traj, litacu_100.atraj, litacu_100.tess)"
    echo "Example: $0 litacu_100.traj (will look for litacu_100.traj, litacu_100.atraj, litacu_100.tess)"
    exit 1
fi

# Get the input filename and extract base name
INPUT_FILENAME="$1"

# Remove extension if present to get base filename
if [[ "$INPUT_FILENAME" == *.traj ]]; then
    BASE_FILENAME="${INPUT_FILENAME%.traj}"
elif [[ "$INPUT_FILENAME" == *.atraj ]]; then
    BASE_FILENAME="${INPUT_FILENAME%.atraj}"
elif [[ "$INPUT_FILENAME" == *.tess ]]; then
    BASE_FILENAME="${INPUT_FILENAME%.tess}"
else
    BASE_FILENAME="$INPUT_FILENAME"
fi

# Construct the full filenames
TRAJ_FILE="${BASE_FILENAME}.traj"
ATRAJ_FILE="${BASE_FILENAME}.atraj"
TESS_FILE="${BASE_FILENAME}.tess"

echo "Running benchmarks for: $BASE_FILENAME"
echo "Files: $TRAJ_FILE, $ATRAJ_FILE, $TESS_FILE"
echo ""

# Read benchmarks
echo "=== Read Benchmarks ==="
echo "Reading .traj file:"
python -m timeit -s "from atomict.io.formats.traj import read_traj" "read_traj('$TRAJ_FILE')"

echo "Reading .atraj file:"
python -m timeit -s "from atomict.io.msgpack import load_msgpack_trajectory" "load_msgpack_trajectory('$ATRAJ_FILE')"

echo "Reading .tess file:"
python -m timeit -s "from atomict.io.msgpack import load_msgpack_trajectory" "load_msgpack_trajectory('$TESS_FILE')"

echo ""
echo "=== Write Benchmarks ==="
echo "Writing .traj file:"
python -m timeit -n 3 -s "from atomict.io.formats.traj import read_traj, write_traj; atoms, metadata = read_traj('$TRAJ_FILE')" "write_traj(atoms, 'test_output.traj', metadata)"

echo "Writing .atraj file:"
python -m timeit -n 3 -s "from atomict.io.formats.traj import read_traj; from atomict.io.msgpack import save_msgpack_trajectory; atoms, metadata = read_traj('$TRAJ_FILE')" "save_msgpack_trajectory(atoms, 'test_output.atraj', metadata)"

echo "Writing .tess file:"
python -m timeit -n 3 -s "from atomict.io.formats.traj import read_traj; from atomict.io.msgpack import save_msgpack_trajectory; atoms, metadata = read_traj('$TRAJ_FILE')" "save_msgpack_trajectory(atoms, 'test_output.tess', metadata)"

