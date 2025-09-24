from atomict.io.msgpack import save_msgpack_trajectory, load_msgpack_trajectory
from datetime import datetime
import os
import time

FILENAME = "equilibration.atraj"


def get_file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024.0 * 1024.0)


def time_call(func, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    return result, elapsed


def ensure_companion_source(base_atraj_path: str, base_tess_path: str) -> None:
    # Ensure we have a tess source to read from when starting with only an atraj input
    if not os.path.exists(base_tess_path):
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
    base_path = os.path.abspath(FILENAME)
    base_dir = os.path.dirname(base_path)
    stem, src_ext = os.path.splitext(os.path.basename(base_path))

    base_atraj = base_path if src_ext == '.atraj' else os.path.join(base_dir, f"{stem}.atraj")
    base_tess = os.path.join(base_dir, f"{stem}.tess")

    if src_ext == '.atraj' and not os.path.exists(base_atraj):
        raise FileNotFoundError(base_atraj)

    # Make sure we have both source types available
    ensure_companion_source(base_atraj, base_tess)

    # Destination filenames for each conversion
    out_atraj_from_atraj = os.path.join(base_dir, f"{stem}.from_atraj.to_atraj.atraj")
    out_tess_from_atraj = os.path.join(base_dir, f"{stem}.from_atraj.to_tess.tess")
    out_atraj_from_tess = os.path.join(base_dir, f"{stem}.from_tess.to_atraj.atraj")
    out_tess_from_tess = os.path.join(base_dir, f"{stem}.from_tess.to_tess.tess")

    results = []
    results.append(benchmark_conversion(base_atraj, out_atraj_from_atraj))
    results.append(benchmark_conversion(base_atraj, out_tess_from_atraj))
    results.append(benchmark_conversion(base_tess, out_atraj_from_tess))
    results.append(benchmark_conversion(base_tess, out_tess_from_tess))

    # Include a brief run timestamp header
    print(f"Benchmark run: {datetime.now().isoformat(timespec='seconds')}")
    print_report(results)


if __name__ == '__main__':
    main()
