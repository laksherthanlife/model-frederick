from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import tomllib

from provision_public_inputs import external_directory, nonsymlink_path, tracked_paths


ROOT = Path(__file__).resolve().parents[1]
PORTABLE_MANIFEST = "data/frozen_evidence/native_v1/manifest.json"
PORTABLE_SHA256 = "eff317f689fe096f6b3df35e050f557b51cd7e51cdb192bb95775eef8ef529ec"
PUBLIC_ENVIRONMENT = {
    "YSTWIN_YEAST_GEM": "data/gem/yeast-GEM.xml.gz",
    "YSTWIN_EC_YEAST_GEM": "data/gem/ecYeastGEM_batch.xml.gz",
    "YSTWIN_THERMO": "data/thermo",
}
SOURCE_TESTS = (
    "tests/test_kinetic_sbml.py::test_published_hog_source_native_values_simple_rate_and_short_integration",
    "tests/test_kinetic_sbml.py::test_published_hog_trajectory_agrees_between_implicit_solvers",
    "tests/test_kinetic_sbml.py::test_published_hog_reaction_rates_match_independent_roadrunner",
    "tests/test_native_reference_models.py",
)
SOURCE_RUNTIME_ARGUMENTS = {
    "octave": ("--no-gui", "--quiet", "--no-history", "--no-init-file", "--no-site-file", "--eval",
               "printf('%s\\n', jsonencode(struct('version', version)));"),
    "pdftotext": ("-v",),
}


@dataclass(frozen=True)
class Check:
    name: str
    argv: tuple[str, ...]


def check_environment(root=ROOT):
    config = tomllib.loads((Path(root) / "pyproject.toml").read_text(encoding="utf-8"))
    quality = config["tool"]["ystwin"]["quality"]
    expected_python = quality["python-version"]
    errors = []
    if platform.python_version() != expected_python:
        errors.append(f"Python must be {expected_python}, not {platform.python_version()}")
    requirements = list(config["project"]["dependencies"])
    for extra in quality["extras"]:
        requirements.extend(config["project"]["optional-dependencies"][extra])
    versions = {}
    for requirement in requirements:
        name, separator, expected = requirement.partition("==")
        if not separator or not expected or any(symbol in expected for symbol in ";<>=, "):
            errors.append(f"quality dependencies must be exactly pinned: {requirement}")
            continue
        try:
            actual = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            actual = None
        versions[name] = {"expected": expected, "actual": actual}
        if actual != expected:
            errors.append(f"{name}: expected {expected}, got {actual}")
    for module in ("numpy", "scipy", "pandas", "cobra", "optlang", "swiglpk", "sympy",
                   "pytfa", "equilibrator_api", "libsbml", "roadrunner"):
        try:
            importlib.import_module(module)
        except Exception as error:
            errors.append(f"required module {module} could not be imported: {error}")
    solver = None
    try:
        import cobra
        import swiglpk
        from optlang import glpk_interface

        from ystwin.fba.solver import configure

        model = cobra.Model("quality_solver_probe")
        settings = configure(model)
        if not isinstance(model.solver, glpk_interface.Model):
            errors.append("the configured numerical solver is not optlang's GLPK interface")
        solver = {"interface": model.solver.interface.__name__, "glpk": swiglpk.glp_version(),
                  "tolerance": settings.tolerance}
    except Exception as error:
        errors.append(f"required GLPK solver could not be configured: {error}")
    return {"passed": not errors, "python": platform.python_version(), "expected_python": expected_python,
            "platform": platform.platform(), "dependencies": versions, "solver": solver, "errors": errors}


def _source_runtime_version(name, completed):
    if name == "octave":
        payload = json.loads(completed.stdout)
        version = payload.get("version") if isinstance(payload, dict) else None
    else:
        # pdftotext exposes its version on stderr, not as JSON.
        banner = completed.stderr.partition("\n")[0]
        prefix = "pdftotext version "
        version = banner.removeprefix(prefix).strip() if banner.startswith(prefix) else None
    if not isinstance(version, str) or not version.strip():
        raise ValueError(f"{name} did not report a version")
    return version


def check_source_runtimes(root=ROOT):
    """Require executable public-source checks, even for a source-only Python adapter.

    Ubuntu CI enforces the pinned upstream versions installed by quality.yml. Local
    non-Linux runs record their actual versions, without claiming Linux reproduction.
    No source files or local Octave startup/history files are read or written here.
    """
    config = tomllib.loads((Path(root) / "pyproject.toml").read_text(encoding="utf-8"))
    pins = config["tool"]["ystwin"]["quality"]["source-runtimes"]
    system = platform.system()
    runtimes, errors = {}, []
    for name, arguments in SOURCE_RUNTIME_ARGUMENTS.items():
        pin = pins[name]
        executable = shutil.which(name)
        expected = pin["linux-version"] if system == "Linux" else None
        row = {"executable": executable, "actual_version": None, "expected_version": expected,
               "ubuntu_package": f"{pin['ubuntu-package']}={pin['ubuntu-package-version']}",
               "returncode": None}
        runtimes[name] = row
        if executable is None:
            row["error"] = f"{name} is required for public source checks; install {row['ubuntu_package']} on Ubuntu 24.04"
        else:
            row["argv"] = [executable, *arguments]
            try:
                completed = subprocess.run(row["argv"], capture_output=True, text=True, check=False, timeout=30)
                row.update(returncode=completed.returncode, stdout=completed.stdout, stderr=completed.stderr)
                if completed.returncode != 0:
                    row["error"] = f"{name} version probe failed with exit {completed.returncode}"
                else:
                    row["actual_version"] = _source_runtime_version(name, completed)
                    if expected is not None and row["actual_version"] != expected:
                        row["error"] = f"{name}: expected {expected}, got {row['actual_version']}"
            except (OSError, subprocess.TimeoutExpired, ValueError) as error:
                row["error"] = f"{name} version probe failed: {error}"
        row["passed"] = "error" not in row
        if not row["passed"]:
            errors.append(row["error"])
    return {"passed": not errors, "platform": system, "runtimes": runtimes, "errors": errors}


def quality_environment(root, output, cache_dir):
    environment = dict(os.environ)
    environment.update({
        "CI": "true",
        "PYTHONPATH": os.pathsep.join((str(root / "src"), str(root))),
        "PYTHONHASHSEED": "0",
        "PYTHONUNBUFFERED": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
        "PYTEST_ADDOPTS": "",
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(output / "matplotlib"),
        "RUFF_CACHE_DIR": str(output / "ruff-cache"),
        "XDG_CACHE_HOME": str(cache_dir / "xdg"),
        "TMPDIR": str(output / "tmp"),
    })
    environment.pop("YSTWIN_OUTPUTS", None)
    environment.update({name: str(root / path) for name, path in PUBLIC_ENVIRONMENT.items()})
    return environment


def _metadata_hashes(root):
    result = {}
    for relative in sorted(tracked_paths(root)):
        if any(part.endswith(".egg-info") for part in Path(relative).parts):
            path = nonsymlink_path(root / relative)
            result[relative] = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return result


def prepare_installation(root, stage):
    root = Path(root).resolve(strict=True)
    stage = external_directory(stage, root)
    stage.mkdir(parents=True, exist_ok=False)
    project = nonsymlink_path(root / "pyproject.toml")
    shutil.copyfile(project, stage / "pyproject.toml")
    source_root = nonsymlink_path(root / "src")
    if not source_root.is_dir():
        raise FileNotFoundError("project installation requires src/")
    copied = ["pyproject.toml"]
    for source in sorted(source_root.rglob("*")):
        relative = source.relative_to(root)
        if (any(part.endswith(".egg-info") or part == "__pycache__" for part in relative.parts)
                or source.suffix == ".pyc"):
            continue
        nonsymlink_path(source)
        if not source.is_file():
            continue
        destination = stage / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        copied.append(relative.as_posix())
    return copied


def install_project(root, output):
    root = Path(root).resolve(strict=True)
    output = external_directory(output, root)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("installation reports require a fresh or empty external directory")
    output.mkdir(parents=True, exist_ok=True)
    before = _metadata_hashes(root)
    stage = output / "source"
    copied = prepare_installation(root, stage)
    config = tomllib.loads((stage / "pyproject.toml").read_text(encoding="utf-8"))
    extras = ",".join(config["tool"]["ystwin"]["quality"]["extras"])
    command = [sys.executable, "-m", "pip", "install", "--disable-pip-version-check", "--no-deps",
               "--no-build-isolation", "--no-index", "--report", str(output / "pip-report.json"),
               f"{stage}[{extras}]"]
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    (output / "tmp").mkdir()
    environment["TMPDIR"] = str(output / "tmp")
    completed = subprocess.run(command, cwd=stage, env=environment, check=False)
    after = _metadata_hashes(root)
    unchanged = before == after
    report = {"schema_version": 1, "passed": completed.returncode == 0 and unchanged,
              "pip_argv": command, "pip_returncode": completed.returncode,
              "source_copy": str(stage), "copied_inputs": copied,
              "omitted_generated_metadata": sorted(before), "metadata_before": before,
              "metadata_after": after, "tracked_metadata_unchanged": unchanged}
    (output / "installation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)
    return completed.returncode if completed.returncode != 0 else (0 if unchanged else 1)


def build_checks(root, output, cache_dir, *, install=False):
    python = sys.executable
    checks = [
        Check("worktree-before", ("git", "diff", "--exit-code", "--")),
        Check("index-before", ("git", "diff", "--cached", "--exit-code", "--")),
    ]
    if install:
        checks.extend([
            Check("install-dependencies", (python, "-m", "pip", "install", "--disable-pip-version-check",
                                           "--require-hashes", "-r", "requirements-quality.txt")),
            Check("install", (python, "scripts/verify_quality.py", "--install-project", "--root", str(root),
                              "--output-dir", str(output / "installation"))),
        ])
    checks.extend([
        Check("dependency-consistency", (python, "-m", "pip", "check")),
        Check("environment", (python, "scripts/verify_quality.py", "--check-environment")),
        Check("source-runtimes", (python, "scripts/verify_quality.py", "--check-source-runtimes")),
        Check("public-inputs", (python, "scripts/provision_public_inputs.py", "--cache-dir", str(cache_dir),
                                "--output-dir", str(output / "provision"))),
        Check("pytest", (python, "-m", "pytest", "-ra", "--require-public-inputs",
                         "-o", f"cache_dir={output / 'pytest-cache'}", "--junitxml", str(output / "junit.xml"))),
        Check("ruff", (python, "-m", "ruff", "check", ".")),
        Check("claims", (python, "scripts/audit_claims.py", "--output-dir", str(output / "claims"))),
        Check("reproducibility", (python, "scripts/audit_reproducibility.py", "--output-dir", str(output / "reproducibility"))),
        Check("portable-native-replay", (python, "scripts/replay_portable_evidence.py", "--historical", "--root", str(root),
                                         "--manifest", PORTABLE_MANIFEST, "--manifest-sha256", PORTABLE_SHA256,
                                         "--output", str(output / "portable-native-replay.json"))),
        Check("source-model-tests", (python, "-m", "pytest", "-ra", *SOURCE_TESTS, "--require-public-inputs",
                                     "-o", f"cache_dir={output / 'source-pytest-cache'}",
                                     "--junitxml", str(output / "source-model-junit.xml"))),
        Check("ecmodel-source", (python, "scripts/ecmodel_isoprenoid_prior.py", "--model",
                                 str(root / PUBLIC_ENVIRONMENT["YSTWIN_EC_YEAST_GEM"]),
                                 "--output-dir", str(output / "ecmodel-source"))),
        Check("stress-map-source", (python, "scripts/audit_stress_map.py", "--output-dir", str(output / "stress-map-source"))),
        Check("hog-source", (python, "scripts/run_hog_learning.py", "--source-dir", str(root / "data/hog2013"),
                             "--output-dir", str(output / "hog-source"), "--gem-check")),
        Check("public-input-integrity", (python, "scripts/provision_public_inputs.py", "--verify-only",
                                         "--output-dir", str(output / "input-integrity"))),
        Check("worktree-after", ("git", "diff", "--exit-code", "--")),
        Check("index-after", ("git", "diff", "--cached", "--exit-code", "--")),
    ])
    return checks


def run_checks(root, output, checks, environment):
    output = external_directory(output, root)
    if not checks or len({check.name for check in checks}) != len(checks):
        raise ValueError("quality checks must be nonempty and uniquely named")
    for check in checks:
        if not check.name or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for c in check.name):
            raise ValueError("invalid quality check name")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("quality reports require a fresh or empty external directory")
    output.mkdir(parents=True, exist_ok=True)
    (output / "tmp").mkdir()
    results = []
    for check in checks:
        started = time.monotonic()
        log_path = output / f"{check.name}.log"
        row = {"name": check.name, "argv": list(check.argv), "returncode": None, "log": log_path.name}
        print(f"[{check.name}] starting: {json.dumps(list(check.argv))}", flush=True)
        with log_path.open("x", encoding="utf-8") as log:
            try:
                completed = subprocess.run(check.argv, cwd=root, env=environment, stdout=log, stderr=subprocess.STDOUT, check=False)
                row["returncode"] = completed.returncode
            except OSError as error:
                row["launch_error"] = str(error)
                log.write(f"Could not start child process: {error}\n")
        row["seconds"] = round(time.monotonic() - started, 3)
        row["passed"] = row["returncode"] == 0
        results.append(row)
        summary = {"schema_version": 1, "passed": all(result["passed"] for result in results),
                   "complete": len(results) == len(checks), "checks": results}
        (output / "status.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"[{check.name}] returncode={row['returncode']}; log={log_path}", flush=True)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run every blocking quality check and preserve each child's actual exit code")
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--install", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-environment", action="store_true")
    mode.add_argument("--check-source-runtimes", action="store_true")
    mode.add_argument("--install-project", action="store_true")
    args = parser.parse_args(argv)
    if args.install_project:
        if args.output_dir is None:
            parser.error("--install-project requires an external --output-dir")
        return install_project(args.root, args.output_dir)
    if args.check_environment or args.check_source_runtimes:
        report = check_environment(args.root) if args.check_environment else check_source_runtimes(args.root)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if report["passed"] else 1
    if args.output_dir is None or args.cache_dir is None:
        parser.error("--output-dir and --cache-dir are required outside the checkout")
    root = args.root.resolve(strict=True)
    output = external_directory(args.output_dir, root)
    cache_dir = external_directory(args.cache_dir, root)
    if output.is_relative_to(cache_dir) or cache_dir.is_relative_to(output):
        parser.error("reports and the reusable public-input cache must be separate")
    checks = build_checks(root, output, cache_dir, install=args.install)
    report = run_checks(root, output, checks, quality_environment(root, output, cache_dir))
    print(f"Quality: {'PASS' if report['passed'] else 'FAIL'}; diagnostics: {output}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
