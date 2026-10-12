"""Install pinned RNA annotators locally; requires Python 3.12 and native compilers."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPOS = {
    "RNAView": ("https://github.com/rcsb/RNAView.git", "master", "6349bbd94584baecf54161f40206308471be665c"),
    "bpnet": ("https://github.com/computational-biology/bpnet.git", "master", "6fd9ebcf8765b004913629ff27aa90143affc428"),
    "fr3d-python": ("https://github.com/BGSU-RNA/fr3d-python.git", "latest", "3e6f92a739ccda126a44969f082dbf6c0eb3ad9f"),
}


def run(command, cwd=None, log=None, env=None):
    command = list(map(str, command))
    print("Running:", " ".join(command), flush=True)
    if log is None:
        subprocess.run(command, cwd=cwd, env=env, check=True)
    else:
        with log.open("a") as handle:
            result = subprocess.run(command, cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Command failed; see {log}")


def locate(names, description):
    for name in names:
        if name and (found := shutil.which(str(name))):
            return str(Path(found).absolute())
    raise RuntimeError(f"Missing {description}; tried {', '.join(str(n) for n in names if n)}")


def compiler(prefix):
    # Apple clang cannot compile this BPNET release's unqualified stack type.
    names = [f"{prefix}-{n}" for n in range(16, 10, -1)] + [prefix]
    for name in names:
        found = shutil.which(name)
        if found:
            version = subprocess.check_output([found, "--version"], text=True)
            if "clang" not in version.lower():
                return found
    raise RuntimeError(f"GNU {prefix} required (macOS: brew install gcc)")


def check(config, env):
    for key in ("python", "parse", "dssr", "rnaview", "bpnet"):
        if not os.access(config[key], os.X_OK):
            raise RuntimeError(f"Missing executable {key}: {config[key]}")
    for key in ("rnaview_home", "bpnet_sysfiles"):
        if not Path(config[key]).is_dir():
            raise RuntimeError(f"Missing reference directory: {config[key]}")
    run([config["python"], "-m", "pip", "check"], env=env)
    run([config["python"], "-c",
         "import barnaba, rnapolis, fr3d, gemmi, mdtraj, numpy; "
         "assert hasattr(numpy, 'chararray'), 'barnaba needs the pinned NumPy release'"], env=env)
    run([config["python"], "-c",
         "import pulp; assert pulp.LpSolverDefault is not None, 'RNApolis CLI requires the pinned CBC solver'"], env=env)
    print("Installation checks passed. Use run_all.sh --pilot to test actual annotations.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--python", help="Python 3.12 executable for the isolated environment")
    ap.add_argument("--dssr", help="Path to your existing licensed x3dna-dssr executable")
    ap.add_argument("--tools-dir", type=Path, default=ROOT / ".runtime-tools")
    ap.add_argument("--parse-source", type=Path, default=ROOT / "PARSE-RNA",
                    help="PARSE engine checkout (default: PARSE-RNA/, cloned as in the main README)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true", help="Verify the existing installation without rebuilding")
    args = ap.parse_args()
    dest = args.tools_dir.resolve()
    env = {**os.environ, "PIP_CACHE_DIR": str(dest / "pip-cache"), "PIP_DISABLE_PIP_VERSION_CHECK": "1"}
    config_path = dest / "config.json"
    if args.check:
        check(json.loads(config_path.read_text()), env)
        return
    if args.jobs < 1:
        ap.error("--jobs must be positive")
    python = locate([args.python] if args.python else ["python3.12", sys.executable], "Python 3.12")
    version = subprocess.check_output([python, "-c", "import sys; print('.'.join(map(str,sys.version_info[:2])))"], text=True).strip()
    if version != "3.12":
        ap.error(f"Use --python /path/to/python3.12; found {version}")
    dssr = locate([args.dssr] if args.dssr else ["x3dna-dssr", Path(python).parent / "x3dna-dssr"], "licensed DSSR binary (--dssr)")
    cc, cxx = compiler("gcc"), compiler("g++")
    fc = locate(["gfortran", "gfortran-15", "gfortran-14", "gfortran-13"], "GNU Fortran")
    for name in ("git", "make", "cmake"):
        locate([name], name)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "src").mkdir(exist_ok=True)
    for name, (url, branch, commit) in REPOS.items():
        checkout = dest / "src" / name
        if not checkout.exists():
            run(["git", "clone", "--branch", branch, url, checkout])
        head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=checkout, text=True).strip()
        dirty = subprocess.check_output(["git", "diff", "HEAD", "--name-only"], cwd=checkout, text=True).strip()
        if dirty:
            raise RuntimeError(f"Tracked source changes in {checkout}; preserve them and use a fresh --tools-dir")
        if head != commit:
            run(["git", "fetch", "origin", commit], cwd=checkout)
            run(["git", "checkout", "--detach", commit], cwd=checkout)
    vpython = dest / "venv/bin/python"
    if not vpython.exists():
        run([python, "-m", "venv", dest / "venv"])
    run([vpython, "-m", "pip", "install", "-r", HERE / "requirements.lock.txt"], env=env)
    run([vpython, "-m", "pip", "install", "--no-deps", dest / "src/fr3d-python"], env=env)
    rv = dest / "src/RNAView"
    run(["make", "-B", f"-j{args.jobs}", f"CC={cc}", "CFLAGS=-I./include -O3"], cwd=rv, log=dest / "rnaview-build.log")
    bp = dest / "src/bpnet/src"
    log = dest / "bpnet-build.log"
    run([fc, "-O3", "-std=legacy", "-fallow-argument-mismatch", "-c", "libbpfind.f", "-o", "libbpfind.o"], cwd=bp, log=log)
    run([cxx, "-O3", "-std=c++11", "-I./include", "-c", "main.cpp", "-o", "main.o"], cwd=bp, log=log)
    run([cxx, "main.o", "libbpfind.o", "-o", "bpnet", "-lgfortran"], cwd=bp, log=log)
    source = args.parse_source.resolve()
    build = dest / "parse-build"
    command = ["cmake", "-S", source, "-B", build, "-DCMAKE_BUILD_TYPE=Release"]
    gemmi = source / "build/_deps/gemmi-src"
    if gemmi.exists():
        command.append(f"-DFETCHCONTENT_SOURCE_DIR_GEMMI={gemmi}")
    run(command, log=dest / "parse-build.log")
    run(["cmake", "--build", build, "--target", "parse", f"-j{args.jobs}"], log=dest / "parse-build.log")
    packages = json.loads(subprocess.check_output([vpython, "-m", "pip", "list", "--format=json"], env=env, text=True))
    config = {"python": str(vpython), "parse": str(build / "parse"), "dssr": dssr,
              "rnaview_home": str(rv), "rnaview": str(rv / "bin/rnaview"),
              "bpnet": str(bp / "bpnet"), "bpnet_sysfiles": str(bp.parent / "sysfiles"),
              "packages": {r["name"]: r["version"] for r in packages},
              "git_commits": {name: values[2] for name, values in REPOS.items()},
              "compilers": {name: subprocess.check_output([exe, "--version"], text=True).splitlines()[0]
                            for name, exe in (("c", cc), ("cxx", cxx), ("fortran", fc))}}
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    check(config, env)
    print(f"Configuration: {config_path}")


if __name__ == "__main__":
    main()
