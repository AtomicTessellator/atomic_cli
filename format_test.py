from atomict.io.msgpack import save_msgpack_trajectory, load_msgpack_trajectory
from atomict.io.formats.traj import read_traj, write_traj
from datetime import datetime
import os
import time
import sys


def get_file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024.0 * 1024.0)


def time_call(func, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    return result, elapsed


def ensure_companion_source(base_atraj_path: str, base_tess_path: str) -> None:
    # Retained for compatibility in case callers rely on generating a tess source
    if not os.path.exists(base_tess_path) and os.path.exists(base_atraj_path):
        atoms, metadata = load_msgpack_trajectory(base_atraj_path)
        save_msgpack_trajectory(atoms, base_tess_path, metadata)


def benchmark_conversion(src_path: str, dst_path: str):
    (atoms, metadata), load_s = time_call(load_msgpack_trajectory, src_path)
    _, save_s = time_call(save_msgpack_trajectory, atoms, dst_path, metadata)
    src_mb = get_file_size_mb(src_path)
    dst_mb = get_file_size_mb(dst_path)
    total_s = load_s + save_s
    return {
        'case': f"{os.path.splitext(src_path)[1][1:]}→{os.path.splitext(dst_path)[1][1:]}",
        'src': os.path.basename(src_path),
        'dst': os.path.basename(dst_path),
        'load_s': load_s,
        'save_s': save_s,
        'total_s': total_s,
        'src_mb': src_mb,
        'dst_mb': dst_mb,
    }

def benchmark_traj_same(src_traj_path: str, dst_traj_path: str):
    (atoms, metadata), load_s = time_call(read_traj, src_traj_path)
    _, save_s = time_call(write_traj, atoms, dst_traj_path, metadata)
    src_mb = get_file_size_mb(src_traj_path)
    dst_mb = get_file_size_mb(dst_traj_path)
    total_s = load_s + save_s
    return {
        'case': 'traj→traj',
        'src': os.path.basename(src_traj_path),
        'dst': os.path.basename(dst_traj_path),
        'load_s': load_s,
        'save_s': save_s,
        'total_s': total_s,
        'src_mb': src_mb,
        'dst_mb': dst_mb,
    }

def benchmark_tess_same(src_tess_path: str, dst_tess_path: str):
    (atoms, metadata), load_s = time_call(load_msgpack_trajectory, src_tess_path)
    _, save_s = time_call(save_msgpack_trajectory, atoms, dst_tess_path, metadata)
    src_mb = get_file_size_mb(src_tess_path)
    dst_mb = get_file_size_mb(dst_tess_path)
    total_s = load_s + save_s
    return {
        'case': 'tess→tess',
        'src': os.path.basename(src_tess_path),
        'dst': os.path.basename(dst_tess_path),
        'load_s': load_s,
        'save_s': save_s,
        'total_s': total_s,
        'src_mb': src_mb,
        'dst_mb': dst_mb,
    }



def print_report(rows):
    headers = (
        ("CASE", 12),
        ("SRC", 28),
        ("DST", 32),
        ("LOAD(s)", 10),
        ("SAVE(s)", 10),
        ("TOTAL(s)", 10),
        ("SRC(MB)", 10),
        ("DST(MB)", 10),
        ("DST/SRC", 9),
    )
    fmt = (
        f"{{case:<{headers[0][1]}.{headers[0][1]}}}"
        f"{{src:<{headers[1][1]}.{headers[1][1]}}}"
        f"{{dst:<{headers[2][1]}.{headers[2][1]}}}"
        f"{{load_s:>{headers[3][1]}.3f}}"
        f"{{save_s:>{headers[4][1]}.3f}}"
        f"{{total_s:>{headers[5][1]}.3f}}"
        f"{{src_mb:>{headers[6][1]}.2f}}"
        f"{{dst_mb:>{headers[7][1]}.2f}}"
        f"{{ratio:>{headers[8][1]}.2f}}"
    )
    line = ''.join(h[:w].ljust(w) for h, w in headers)
    sep = '-' * len(line)
    print(sep)
    print(line)
    print(sep)
    for r in rows:
        r['ratio'] = (r['dst_mb'] / r['src_mb']) if r['src_mb'] > 0 else 0.0
        print(fmt.format(**r))
    print(sep)


def main():
    # Accept optional base name on CLI: e.g. `python format_test.py 800K`
    arg = sys.argv[1] if len(sys.argv) > 1 else FILENAME

    base_path = os.path.abspath(arg)
    base_dir = os.path.dirname(base_path)
    stem, ext = os.path.splitext(os.path.basename(base_path))

    # If user provided a stem without extension, keep it; otherwise remove extension from stem
    if ext:
        base_stem = stem
    else:
        base_stem = os.path.basename(base_path)

    traj_path = os.path.join(base_dir, f"{base_stem}.traj")
    tess_path = os.path.join(base_dir, f"{base_stem}.tess")

    if not os.path.exists(traj_path):
        raise FileNotFoundError(traj_path)
    if not os.path.exists(tess_path):
        raise FileNotFoundError(tess_path)

    out_traj_from_traj = os.path.join(base_dir, f"{base_stem}.from_traj.to_traj.traj")
    out_tess_from_tess = os.path.join(base_dir, f"{base_stem}.from_tess.to_tess.tess")

    results = []
    # results.append(benchmark_traj_same(traj_path, out_traj_from_traj))
    results.append(benchmark_tess_same(tess_path, out_tess_from_tess))

    print(f"Benchmark run: {datetime.now().isoformat(timespec='seconds')}")
    print_report(results)


if __name__ == '__main__':
    main()
