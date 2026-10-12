"""Sequential, resumable runtime benchmark of all seven RNA annotators.

Primary time is fresh-worker wall time, including its common Python launcher,
tool startup/imports, required conversions, annotation, and output validation.
No network access or cached annotations are used by the adapters.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import shutil
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TOOLS = ("parse", "parse_noscore", "dssr", "fr3d", "rnaview", "barnaba", "rnapolis", "bpnet")
THREAD_ENV = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS")}


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, data):
    temp = path.with_suffix(path.suffix + ".part")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def choose_structures(rows, limit, seed):
    """Select before timing, balanced across four predeclared size strata."""
    if not limit or limit >= len(rows):
        return sorted(rows, key=lambda r: r["pdb"])
    if limit < 1:
        raise ValueError("limit must be nonnegative")
    buckets = [[] for _ in range(4)]
    for row in rows:
        size = row["nucleotides"]
        buckets[0 if size <= 50 else 1 if size <= 200 else 2 if size <= 1000 else 3].append(row)
    rng = random.Random(seed)
    for bucket in buckets:
        rng.shuffle(bucket)
    chosen = []
    while len(chosen) < limit:
        for bucket in buckets:
            if bucket and len(chosen) < limit:
                chosen.append(bucket.pop())
    return sorted(chosen, key=lambda r: r["pdb"])


def prepare_input(source, target):
    """Keep the first deposited model, preserving mmCIF categories and coordinates."""
    import gemmi
    doc = gemmi.cif.read_file(str(source))
    block = doc.sole_block()
    table = block.find_mmcif_category("_atom_site.")
    if not len(table):
        raise ValueError(f"No atom records: {source}")
    loop = table.loop
    model_tag = "_atom_site.pdbx_PDB_model_num"
    if model_tag in loop.tags:
        col = list(loop.tags).index(model_tag)
        first_model = table[0][col]
        rows = [list(row) for row in table if row[col] == first_model]
        if len(rows) != len(table):
            loop.set_all_values(list(map(list, zip(*rows))))
            doc.write_file(str(target))
        else:
            shutil.copyfile(source, target)
    else:
        shutil.copyfile(source, target)


def collect_manifest(args, output):
    with args.cohort.open() as handle:
        rows = [{"pdb": r["PDB"].upper(), "nucleotides": int(r["Nucleotides"])}
                for r in csv.DictReader(handle)]
    if len({r["pdb"] for r in rows}) != len(rows):
        raise ValueError("Duplicate PDB IDs in cohort")
    if args.pdbs:
        wanted = set(args.pdbs.upper().split(","))
        rows = [r for r in rows if r["pdb"] in wanted]
        if {r["pdb"] for r in rows} != wanted:
            raise ValueError("Some --pdbs are absent from the cohort")
    rows = choose_structures(rows, args.limit, args.seed)
    input_dir = output / "inputs"
    input_dir.mkdir(exist_ok=True)
    for row in rows:
        source = args.input / f"{row['pdb']}.cif"
        target = input_dir / source.name
        row["source_sha256"] = sha(source)
        prepare_input(source, target)
        row["input_sha256"] = sha(target)
        row["source_path"] = str(source.resolve())
    return rows


def terminate_tree(proc):
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.communicate()


def run_trial(tool, row, repeat, args, output, config):
    work = output / "work" / f"{row['pdb']}_{tool}_{repeat}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    (work / "scratch").mkdir()
    result_file = work / "result.json"
    command = [config["python"], str(HERE / "run_one.py"), "--tool", tool,
               "--input", str(output / "inputs" / f"{row['pdb']}.cif"),
               "--config", str(args.config), "--result", str(result_file)]
    env = {**os.environ, **THREAD_ENV, "TMPDIR": str(work / "scratch"),
           "XDG_CACHE_HOME": str(work / "cache"),
           "MPLCONFIGDIR": str(output / "mpl-cache"), "PYTHONHASHSEED": "0"}
    start = time.perf_counter()
    proc = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, start_new_session=True)
    try:
        log, _ = proc.communicate(timeout=args.timeout)
        elapsed = time.perf_counter() - start
        if proc.returncode == 0 and result_file.is_file():
            result = json.loads(result_file.read_text())
        else:
            result = {"status": "worker_error", "error": f"exit {proc.returncode}"}
    except subprocess.TimeoutExpired as error:
        terminate_tree(proc)
        elapsed = time.perf_counter() - start
        log = error.output or b""
        result = {"status": "timeout"}
    except BaseException:
        terminate_tree(proc)
        raise
    log_path = output / "logs" / f"{row['pdb']}_{tool}_{repeat}.log"
    log_path.write_bytes(log)
    result.update(pdb=row["pdb"], nucleotides=row["nucleotides"], tool=tool,
                  repeat=repeat, warmup=repeat < 0, wall_seconds=elapsed,
                  worker_returncode=proc.returncode, command=command,
                  log=str(log_path.relative_to(output)))
    shutil.rmtree(work)
    return result


def read_trials(path):
    if not path.exists():
        return []
    # A power loss may leave a partial final append. Preserve it for inspection
    # and rerun that trial; malformed complete records must still fail loudly.
    raw = path.read_bytes()
    if raw and not raw.endswith(b"\n"):
        end = raw.rfind(b"\n") + 1
        path.with_suffix(".interrupted-tail").write_bytes(raw[end:])
        path.write_bytes(raw[:end])
        raw = raw[:end]
    return [json.loads(line) for line in raw.splitlines() if line.strip()]


def hardware():
    info = {"cpu_count": os.cpu_count()}
    if sys.platform == "darwin":
        for name, key in (("cpu_model", "machdep.cpu.brand_string"), ("memory_bytes", "hw.memsize")):
            result = subprocess.run(["sysctl", "-n", key], capture_output=True, text=True)
            info[name] = result.stdout.strip() if result.returncode == 0 else "unavailable"
    elif Path("/proc/cpuinfo").exists():
        info["cpu_model"] = next((s.split(":", 1)[1].strip() for s in
            Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")), platform.machine())
        info["memory"] = Path("/proc/meminfo").read_text().splitlines()[0]
    return info


def summarize(trials, manifest, tools, repeats):
    grouped = defaultdict(list)
    for trial in trials:
        if not trial["warmup"]:
            grouped[(trial["pdb"], trial["tool"])].append(trial)
    summaries = []
    for row in manifest:
        for tool in tools:
            runs = grouped[(row["pdb"], tool)]
            good = [r for r in runs if r["status"] == "ok"]
            complete = len(runs) == repeats and len(good) == repeats
            signatures = {r.get("output_sha256") for r in good}
            stable = len(signatures) == 1
            eligible = complete and stable
            summaries.append({"pdb": row["pdb"], "nucleotides": row["nucleotides"], "tool": tool,
                              "attempted": len(runs), "successful": len(good),
                              "complete": complete, "stable_output": stable,
                              "median_seconds": statistics.median(r["wall_seconds"] for r in good) if eligible else "",
                              "min_seconds": min(r["wall_seconds"] for r in good) if eligible else "",
                              "max_seconds": max(r["wall_seconds"] for r in good) if eligible else "",
                              "pair_count": good[0].get("pair_count", "") if eligible else "",
                              "statuses": ";".join(sorted({r["status"] for r in runs}))})
    return summaries


def write_csv(path, rows, fields=None):
    if not rows and fields is None:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def report(output, trials, manifest, tools, repeats):
    summary = summarize(trials, manifest, tools, repeats)
    write_csv(output / "per_structure.csv", summary)
    lookup = {(r["pdb"], r["tool"]): r for r in summary}
    paired = []
    for baseline in ("parse", "parse_noscore"):
        if baseline not in tools:
            continue
        for row in summary:
            ref = lookup[(row["pdb"], baseline)]
            if row["tool"] != baseline and row["median_seconds"] != "" and ref["median_seconds"] != "":
                paired.append({"pdb": row["pdb"], "nucleotides": row["nucleotides"],
                               "tool": row["tool"], "baseline": baseline,
                               "tool_seconds": row["median_seconds"], "parse_seconds": ref["median_seconds"],
                               "speedup": row["median_seconds"] / ref["median_seconds"]})
    write_csv(output / "speedups.csv", paired,
              ["pdb", "nucleotides", "tool", "baseline", "tool_seconds", "parse_seconds", "speedup"])
    coverage = []
    for tool in tools:
        rows = [r for r in summary if r["tool"] == tool]
        speedups = [r["speedup"] for r in paired if r["tool"] == tool and r["baseline"] == "parse"]
        coverage.append({"tool": tool, "structures": len(rows),
                         "complete_stable": sum(r["median_seconds"] != "" for r in rows),
                         "paired_with_parse": len(speedups),
                         "median_speedup_vs_parse": statistics.median(speedups) if speedups else ""})
    write_csv(output / "coverage.csv", coverage)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=ROOT / ".runtime-tools/config.json")
    ap.add_argument("--input", type=Path, default=ROOT / "data" / "structures")
    ap.add_argument("--cohort", type=Path, default=ROOT / "data/benchmark/dssr_PARSE_time_benchmark.csv")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--tools", default=",".join(TOOLS))
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--warmups", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=600)
    ap.add_argument("--limit", type=int, default=0, help="Size-stratified subset; 0 = full cohort")
    ap.add_argument("--pdbs", help="Comma-separated IDs from the cohort")
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    if args.runs < 1 or args.warmups < 0 or args.limit < 0 or args.timeout <= 0:
        ap.error("Invalid run count, warmups, limit, or timeout")
    args.config = args.config.resolve()
    args.input = args.input.resolve()
    tools = args.tools.split(",")
    if len(set(tools)) != len(tools) or not set(tools) <= set(TOOLS):
        ap.error("Unknown or duplicate tool")
    config = json.loads(args.config.read_text())
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()) and not args.resume:
        ap.error("Output already exists; use --resume or choose a new output directory")
    output.mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir(exist_ok=True)
    metadata = {"schema": 1, "tools": tools, "runs": args.runs, "warmups": args.warmups,
                "timeout": args.timeout, "seed": args.seed, "limit": args.limit, "pdbs": args.pdbs,
                "platform": platform.platform(), "machine": platform.machine(),
                "host": platform.node(), "hardware": hardware(),
                "thread_environment": THREAD_ENV, "config": config,
                "cohort_sha256": sha(args.cohort),
                "script_hashes": {name: sha(HERE / name) for name in
                    ("benchmark_all.py", "run_one.py", "_annotators.py")},
                "executable_hashes": {k: sha(config[k]) for k in sorted({"python"} | {
                    "parse" if t.startswith("parse") else t
                    for t in tools if t.startswith("parse") or t in {"dssr", "rnaview", "bpnet"}})},
                "installed_packages": subprocess.check_output(
                    [config["python"], "-m", "pip", "freeze", "--disable-pip-version-check"],
                    env={**os.environ, "PIP_CACHE_DIR": str(ROOT / ".runtime-tools/pip-cache")},
                    text=True).splitlines()}
    meta_path = output / "metadata.json"
    if args.resume and meta_path.exists():
        if json.loads(meta_path.read_text()) != metadata:
            ap.error("Configuration, code, host, or executables changed; use a new output directory")
        manifest = json.loads((output / "manifest.json").read_text())
        for row in manifest:
            if sha(args.input / f"{row['pdb']}.cif") != row["source_sha256"]:
                ap.error(f"Input changed: {row['pdb']}")
            if sha(output / "inputs" / f"{row['pdb']}.cif") != row["input_sha256"]:
                ap.error(f"Prepared input changed: {row['pdb']}")
    else:
        manifest = collect_manifest(args, output)
        write_json(output / "manifest.json", manifest)
        write_json(meta_path, metadata)
    trials_path = output / "trials.jsonl"
    trials = read_trials(trials_path)
    done = {(r["pdb"], r["tool"], r["repeat"]) for r in trials}
    if len(done) != len(trials):
        ap.error("Duplicate trial records")
    print(f"{len(manifest)} structures; {len(tools)} modes; {args.warmups} warm-up + {args.runs} repeats", flush=True)
    started = time.perf_counter()
    try:
        with trials_path.open("a") as handle:
            for index, row in enumerate(manifest, 1):
                rng = random.Random(f"{args.seed}:{row['pdb']}")
                for repeat in range(-args.warmups, args.runs):
                    order = list(tools)
                    rng.shuffle(order)
                    for tool in order:
                        if (row["pdb"], tool, repeat) in done:
                            continue
                        result = run_trial(tool, row, repeat, args, output, config)
                        handle.write(json.dumps(result, sort_keys=True) + "\n")
                        handle.flush()
                        os.fsync(handle.fileno())
                        trials.append(result)
                        print(f"[{index}/{len(manifest)}] {row['pdb']} {tool} run={repeat} "
                              f"{result['status']} {result['wall_seconds']:.3f}s", flush=True)
                report(output, trials, manifest, tools, args.runs)
    finally:
        report(output, trials, manifest, tools, args.runs)
    print(f"Finished in {time.perf_counter()-started:.1f}s: {output}", flush=True)


if __name__ == "__main__":
    main()
