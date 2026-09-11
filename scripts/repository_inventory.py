from __future__ import annotations

import argparse
import ast
from collections import defaultdict
from copy import deepcopy
import csv
import fnmatch
from graphlib import CycleError, TopologicalSorter
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from typing import Any


SCHEMA_VERSION = 1
OWNED_PATHS = (
    "scripts/repository_inventory.py",
    "data/repository_inventory.json",
    "data/implementation_workflow.json",
    "data/deletion_manifest.json",
    "tests/test_repository_inventory.py",
)
DOCUMENTS = (
    "repository_inventory.json",
    "implementation_workflow.json",
    "deletion_manifest.json",
)
ROLES = (
    "active_runtime",
    "supported_public_library",
    "validation_reproduction",
    "recorded_negative_result",
    "superseded",
    "unowned_requires_decision",
)
FROZEN_MANIFEST = "data/frozen_evidence/native_v1/manifest.json"
FROZEN_MANIFEST_SHA256 = "eff317f689fe096f6b3df35e050f557b51cd7e51cdb192bb95775eef8ef529ec"
STATUS_ONLY = frozenset({"S.csv", "jws.html", ".ystwin.env", ".env", ".env.local"})
MECHANICAL_REFERENCES = frozenset(OWNED_PATHS)
CLEANUP_APPROVAL_ID = "reviewed-cleanup-20260908-v1"
CLEANUP_APPROVAL = {
    "actor": "user",
    "date": "2026-09-08",
    "selection": "Approve all nine actions",
    "scope": "Only the nine exact listed actions, after their replacement and historical-integrity preconditions pass.",
}
CLEANUP_APPROVAL_2_ID = "reviewed-cleanup-20260910-partial-orders"
CLEANUP_APPROVAL_2 = {
    "actor": "user",
    "date": "2026-09-10",
    "selection": "Approve removing both superseded partial-order originals",
    "scope": "Only the two exact listed actions. partial_order_v2 is a verified strict superset.",
}
CLEANUP_APPROVAL_3_ID = "reviewed-cleanup-20260911-autofluorescence-supersession"
CLEANUP_APPROVAL_3 = {
    "actor": "user",
    "date": "2026-09-11",
    "selection": "Approve relocating the superseded autofluorescence first pass",
    "scope": "One exact relocation. The content is retained under docs/superseded/, not deleted.",
}
#: ``{approval id: the recorded decision}``. Each contract names the batch that authorised it, so a
#: later approval can never be backdated into an earlier batch's sealed record.
CLEANUP_APPROVALS = {CLEANUP_APPROVAL_ID: CLEANUP_APPROVAL, CLEANUP_APPROVAL_2_ID: CLEANUP_APPROVAL_2,
                     CLEANUP_APPROVAL_3_ID: CLEANUP_APPROVAL_3}
CLEANUP_BASELINE = "84561f51502b5b4bf66c052dc88c853038917934"
RELEASE_SOURCE = "data/native_law_v2/granados/development/phase_c_01/release_phase_c.py"
RELEASE_ASSET = RELEASE_SOURCE + ".source"
TASK_TABLE = "outputs/regulation_task_priority.csv"
TASK_REPORT = "outputs/regulation_matched_tasks.json"
CLEANUP_CONTRACTS = {
    name: {
        "action": action,
        "original": {
            "sha256": sha256, "size_bytes": size, "git_object_id": object_id,
            "git_commit": commit, "tracked": object_id is not None, "ignored": name == "jws.html",
        },
        "replacement_paths": list(replacements),
        "batch": batch,
    }
    for name, action, sha256, size, object_id, commit, replacements, batch in (
        ("src/ystwin.egg-info/PKG-INFO", "remove_generated_file",
         "0f1fff8b361edaf684426c3b2941988b0b7bd3e642cac0c919da536f91fa6c73", 390,
         "b688cea1fd5651dfba7213f20cc368b8e5c1b386", CLEANUP_BASELINE, ("pyproject.toml",), CLEANUP_APPROVAL_ID),
        ("src/ystwin.egg-info/SOURCES.txt", "remove_generated_file",
         "cf81cc0b63e3ecf50bdad3cd94a02e3c1b20759fa92193f018b040c7c1cd7a14", 386,
         "6fe3efee288ff9e16c6b602233704630d23ca30d", CLEANUP_BASELINE,
         ("pyproject.toml", "scripts/repository_inventory.py", "data/repository_inventory.json"), CLEANUP_APPROVAL_ID),
        ("src/ystwin.egg-info/dependency_links.txt", "remove_generated_file",
         "01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b", 1,
         "8b137891791fe96927ad78e64b0aad7bded08bdc", CLEANUP_BASELINE, ("pyproject.toml",), CLEANUP_APPROVAL_ID),
        ("src/ystwin.egg-info/requires.txt", "remove_generated_file",
         "b9a0b9fff15425d30dfb5493736e001854c6be2481d670389168ca8886524997", 94,
         "a98dc2ea5acf4707172e5961016f57c34b54b6c9", CLEANUP_BASELINE,
         ("pyproject.toml", "requirements-quality.txt"), CLEANUP_APPROVAL_ID),
        ("src/ystwin.egg-info/top_level.txt", "remove_generated_file",
         "ca042f34af4a425f49dd6ff3e5f5a313c830ecb52297db7aac5dad79a8fd516a", 7,
         "08019e9d5e06a09e46792af8793e566c4a20973f", CLEANUP_BASELINE,
         ("pyproject.toml", "src/ystwin/__init__.py"), CLEANUP_APPROVAL_ID),
        ("S.csv", "remove_superseded_scratch",
         "b52e7e2388a9d371b4cd40db5bcbddae7de2f6373a79a4d4844b5aea43d4a73d", 2544,
         "316f94b15acd304d9bd4eb9fdc5f72d7e1ff8204", "385b6b8dcf0f4f1a499df3a75de509237c27188f",
         ("src/ystwin/pathway/gem_environment.py", "scripts/gem_environment_predictions.py",
          "outputs/gem_environment_predictions.csv", "outputs/gem_environment_scores.csv", "tests/test_gem_environment.py"), CLEANUP_APPROVAL_ID),
        ("jws.html", "remove_local_download",
         "95827d2849c6aea0c2ed1881218c38b336ca8809686836636ab694a3e0f324f1", 60203, None, None,
         ("data/hog2013/model_wt.xml", "data/hog2013/sources.json",
          "data/native_reference_models/jalihal2021/nutrient-signaling-sbml.xml"), CLEANUP_APPROVAL_ID),
        (TASK_TABLE, "replace_superseded_result",
         "6ed033a695c8d723a67034aa20ec35899c34ec9655772dceadbe9a4b27041ded", 172,
         "764ef3aa9eb9c40f0b642243098d15002874c41d", CLEANUP_BASELINE,
         (TASK_REPORT, "data/artifact_registry.json", "data/current_claims.json",
          "scripts/run_regulation.py", "tests/test_run_regulation.py"), CLEANUP_APPROVAL_ID),
        (RELEASE_SOURCE, "relocate_immutable_source_evidence",
         "68c31da88770e4718bc5668791edd38de6a054fa821aaa074a2ef48d6c676d8a", 42204,
         "fd9dbe4b7afbeec94dd543754396d444b3c45880", CLEANUP_BASELINE, (RELEASE_ASSET,), CLEANUP_APPROVAL_ID),
        ("outputs/partial_order_synthetic_summary.csv", "remove_superseded_result",
         "85f0e481b3f9061ef9eed6ab087ba00a629e4b11d8a64293225a5e43cff05c26", 2758,
         "411c60baeb2f8594a45f0b41f8c9e88fd2fcc873", CLEANUP_BASELINE,
         ("outputs/partial_order_v2/partial_order_synthetic_summary.csv",), CLEANUP_APPROVAL_2_ID),
        ("outputs/partial_order_synthetic_demo.json", "remove_superseded_result",
         "459b70cf347f0354981a69e90033e263d18787f7d606484488b61e1b94676669", 302167,
         "b09fbd94d1d5e985a77d50e84561154442d4e325", CLEANUP_BASELINE,
         ("outputs/partial_order_v2/partial_order_synthetic_demo.json",), CLEANUP_APPROVAL_2_ID),
        # A relocation, not a removal: the content is retained under docs/superseded/ with a header
        # naming what overturned it. The old path leaves the tree, which is what needs approving.
        ("docs/research/AUTOFLUORESCENCE_FOUND.md", "relocate_superseded_document",
         "168b5070451ef55bb3ade2322a9cc32e7564a7056d42aeb99b369ee980047bd0", 4244,
         "bbfd7b1128521f03f340c55baed39865e33c9198", CLEANUP_BASELINE,
         ("docs/superseded/autofluorescence-first-pass.md",), CLEANUP_APPROVAL_3_ID),
    )
}


class InventoryError(ValueError):
    pass


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def serialized(value: Any, level=0) -> str:
    if isinstance(value, dict):
        if not value:
            return "{}"
        if all(not isinstance(item, (dict, list)) for item in value.values()):
            compact = json.dumps(value, allow_nan=False)
            if len(compact) + level <= 120:
                return compact
        pad = " " * (level + 2)
        rows = [pad + json.dumps(key) + ": " + serialized(item, level + 2) for key, item in value.items()]
        return "{\n" + ",\n".join(rows) + "\n" + " " * level + "}"
    if isinstance(value, list):
        if not value:
            return "[]"
        if all(not isinstance(item, (dict, list)) for item in value):
            compact = json.dumps(value, allow_nan=False)
            if len(compact) + level <= 120:
                return compact
        return "[\n" + ",\n".join(" " * (level + 2) + serialized(item, level + 2) for item in value) + "\n" + " " * level + "]"
    return json.dumps(value, allow_nan=False)


def relative_path(value: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or "\0" in value:
        raise InventoryError(f"Invalid repository-relative path: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value or value == ".":
        raise InventoryError(f"Unsafe repository-relative path: {value!r}")
    return value


def safe_file(root: Path, name: str) -> Path:
    path = root / relative_path(name)
    if path.resolve() != path or not path.is_file():
        raise InventoryError(f"Missing or symlinked input: {name}")
    return path


def read_public(root: Path, name: str) -> bytes:
    relative_path(name)
    if name in STATUS_ONLY or any(part.startswith(".env") for part in PurePosixPath(name).parts):
        raise InventoryError(f"Status-only input must not be read: {name}")
    return safe_file(root, name).read_bytes()


def _unique_json(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InventoryError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def parse_json(data: bytes) -> dict:
    try:
        value = json.loads(data, object_pairs_hook=_unique_json)
        canonical(value)
    except (ValueError, TypeError) as exc:
        raise InventoryError(f"Invalid finite JSON: {exc}") from exc
    if not isinstance(value, dict):
        raise InventoryError("Expected a JSON object")
    return value


def git(root: Path, *arguments: str, allowed=(0,)) -> bytes:
    if not arguments or arguments[0] not in {"ls-files", "ls-tree", "rev-parse", "check-ignore", "show"}:
        raise InventoryError("Inventory Git access is read-only; cleanup commands are never executed")
    result = subprocess.run(
        ["git", "--no-replace-objects", "--no-optional-locks", "-C", str(root), *arguments],
        capture_output=True,
        check=False,
    )
    if result.returncode not in allowed:
        raise InventoryError(f"git {' '.join(arguments[:2])} failed: {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def git_catalog(root: Path, baseline: str) -> tuple[dict, dict, dict]:
    if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", baseline):
        raise InventoryError("Baseline must be an exact Git commit object ID")
    current = {}
    for record in git(root, "ls-files", "--stage", "-z").split(b"\0"):
        if not record:
            continue
        metadata, raw_path = record.split(b"\t", 1)
        mode, oid, stage = metadata.decode().split()
        name = relative_path(raw_path.decode("utf-8"))
        if stage != "0":
            raise InventoryError(f"Unmerged Git index path: {name}")
        current[name] = {"mode": mode, "object_id": oid}
    original = {}
    for record in git(root, "ls-tree", "-r", "-z", "--full-tree", baseline).split(b"\0"):
        if not record:
            continue
        metadata, raw_path = record.split(b"\t", 1)
        mode, kind, oid = metadata.decode().split()
        name = relative_path(raw_path.decode("utf-8"))
        original[name] = {"mode": mode, "object_id": oid, "kind": kind}
    metadata = {
        "baseline_commit": baseline,
        "baseline_tree_object_id": git(root, "rev-parse", f"{baseline}^{{tree}}").decode().strip(),
        "observed_head": git(root, "rev-parse", "HEAD").decode().strip(),
        "tracked_pathset_sha256": digest(sorted(current)),
        "baseline_blob_catalog_sha256": digest(original),
        "git_object_hash_policy": "Baseline/index object IDs are Git metadata, not hashes asserted for mutable worktree code.",
    }
    return current, original, metadata


def path_status(root: Path, name: str, tracked: set[str]) -> dict:
    path = root / relative_path(name)
    ignored = bool(git(root, "check-ignore", "--no-index", "--", name, allowed=(0, 1)))
    return {
        "path": name,
        "tracked": name in tracked,
        "exists": path.exists() or path.is_symlink(),
        "symlink": path.is_symlink(),
        "ignored": ignored,
        "status": "tracked" if name in tracked else "untracked_ignored" if ignored else "untracked",
        "contents_read": False,
    }


def module_name(name: str) -> str | None:
    path = PurePosixPath(name)
    if path.suffix != ".py":
        return None
    parts = list(path.with_suffix("").parts)
    if parts[0] == "src":
        parts = parts[1:]
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def context_kind(name: str) -> str:
    if name.startswith("src/ystwin/"):
        return "library"
    if name.startswith("scripts/"):
        return "script"
    if name.startswith("tests/"):
        return "test"
    return "historical_evidence"


def python_allowed(name: str) -> bool:
    return name.endswith(".py") and name.startswith((
        "src/ystwin/", "scripts/", "tests/", "archive/", "data/native_law_v2/", "outputs/",
    ))


def _dotted(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        base = _dotted(node.value)
        return f"{base}.{node.attr}" if base else node.attr
    return ""


def static_path(node: ast.AST, assignments: dict, source: str, seen=frozenset()) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name):
        if node.id == "__file__":
            return source
        if node.id not in seen and node.id in assignments:
            return static_path(assignments[node.id], assignments, source, seen | {node.id})
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Div, ast.Add)):
        left = static_path(node.left, assignments, source, seen)
        right = static_path(node.right, assignments, source, seen)
        if left is not None and right is not None:
            return f"{left}/{right}".lstrip("/") if isinstance(node.op, ast.Div) else left + right
    if isinstance(node, ast.JoinedStr):
        parts = []
        for item in node.values:
            value = static_path(item.value, assignments, source, seen) if isinstance(item, ast.FormattedValue) else static_path(item, assignments, source, seen)
            parts.append(value if value is not None else "*")
        return "".join(parts)
    if isinstance(node, ast.Attribute) and node.attr == "parent":
        value = static_path(node.value, assignments, source, seen)
        if value is not None:
            parent = str(PurePosixPath(value).parent)
            return "" if parent == "." else parent
    if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Attribute) and node.value.attr == "parents":
        value = static_path(node.value.value, assignments, source, seen)
        if value is not None and isinstance(node.slice, ast.Constant) and type(node.slice.value) is int:
            parents = PurePosixPath(value).parents
            if 0 <= node.slice.value < len(parents):
                parent = str(parents[node.slice.value])
                return "" if parent == "." else parent
    if isinstance(node, ast.Call):
        function = _dotted(node.func)
        if function in {"paths.outputs_dir", "outputs_dir"}:
            return "outputs"
        if function in {"paths.data_dir", "data_dir"}:
            return "data"
        if function in {"Path", "pathlib.Path", "str"} and len(node.args) == 1:
            return static_path(node.args[0], assignments, source, seen)
        if isinstance(node.func, ast.Attribute) and node.func.attr in {"resolve", "absolute"}:
            return static_path(node.func.value, assignments, source, seen)
    return None


def assignments_in(tree: ast.AST) -> dict:
    assignments = {}
    ambiguous = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    if target.id in assignments and ast.dump(assignments[target.id]) != ast.dump(node.value):
                        ambiguous.add(target.id)
                    assignments[target.id] = node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None:
            assignments[node.target.id] = node.value
    return {key: value for key, value in assignments.items() if key not in ambiguous}


def analyze_python(sources: dict[str, str], all_paths: set[str]) -> dict:
    modules = {module_name(name): name for name in sources}
    script_aliases = {PurePosixPath(name).stem: name for name in sources if PurePosixPath(name).parent == PurePosixPath("scripts")}
    facts, edges, references, dynamic_gaps, writers = {}, [], [], [], []

    for source, text in sorted(sources.items()):
        try:
            tree = ast.parse(text, filename=source)
        except SyntaxError as exc:
            raise InventoryError(f"Cannot inventory invalid Python {source}:{exc.lineno}: {exc.msg}") from exc
        assignments = assignments_in(tree)
        loads = defaultdict(list)
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
                loads[node.id].append(node.lineno)
        exports = []
        for statement in tree.body:
            if isinstance(statement, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "__all__" for target in statement.targets):
                try:
                    declared = ast.literal_eval(statement.value)
                    if isinstance(declared, (list, tuple)) and all(isinstance(item, str) for item in declared):
                        exports = list(declared)
                except (ValueError, TypeError):
                    pass
        facts[source] = {
            "module": module_name(source),
            "context": context_kind(source),
            "declared_exports": exports,
            "definitions": [{"name": node.name, "line": node.lineno, "kind": "class" if isinstance(node, ast.ClassDef) else "function"} for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))],
            "test_definition_count": sum(isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_") for node in ast.walk(tree)),
            "external_imports": [],
        }
        package = module_name(source) or ""
        if not source.endswith("/__init__.py"):
            package = package.rpartition(".")[0]

        def resolve(dotted: str) -> str | None:
            if dotted in modules:
                return modules[dotted]
            if dotted in script_aliases:
                return script_aliases[dotted]
            sibling = PurePosixPath(source).parent / (dotted.replace(".", "/") + ".py")
            return str(sibling) if str(sibling) in sources else None

        def add_edge(target, node, kind, names, deferred, type_only):
            if target == source:
                return
            edges.append({
                "source": source,
                "target": target,
                "kind": kind,
                "line": node.lineno,
                "bound_names": sorted(names),
                "name_use_lines": sorted({line for name in names for line in loads[name]}),
                "deferred": deferred,
                "type_checking_only": type_only,
            })
            target_module = module_name(target)
            if target_module and target.startswith("src/ystwin/"):
                parts = target_module.split(".")
                for length in range(1, len(parts)):
                    parent = modules.get(".".join(parts[:length]))
                    if parent and parent != source:
                        edges.append({
                            "source": source,
                            "target": parent,
                            "kind": "package_initialization",
                            "line": node.lineno,
                            "bound_names": [],
                            "name_use_lines": [],
                            "deferred": deferred,
                            "type_checking_only": type_only,
                        })

        class Visitor(ast.NodeVisitor):
            def __init__(self):
                self.function_depth = 0
                self.type_only = False

            def visit_FunctionDef(self, node):
                self.function_depth += 1
                self.generic_visit(node)
                self.function_depth -= 1

            visit_AsyncFunctionDef = visit_FunctionDef

            def visit_If(self, node):
                previous = self.type_only
                self.visit(node.test)
                self.type_only = previous or _dotted(node.test) in {"TYPE_CHECKING", "typing.TYPE_CHECKING"}
                for statement in node.body:
                    self.visit(statement)
                self.type_only = previous
                for statement in node.orelse:
                    self.visit(statement)

            def visit_Import(self, node):
                for alias in node.names:
                    target = resolve(alias.name)
                    if target:
                        add_edge(target, node, "import", [alias.asname or alias.name.split(".")[0]], self.function_depth > 0, self.type_only)
                    else:
                        facts[source]["external_imports"].append(alias.name)
                        if alias.name == "ystwin" or alias.name.startswith("ystwin."):
                            dynamic_gaps.append({"path": source, "line": node.lineno, "kind": "unresolved_local_import", "module": alias.name})

            def visit_ImportFrom(self, node):
                base = node.module or ""
                if node.level:
                    parts = package.split(".") if package else []
                    if node.level > len(parts):
                        dynamic_gaps.append({"path": source, "line": node.lineno, "kind": "relative_import_outside_package"})
                        return
                    base = ".".join(parts[:len(parts) - node.level + 1] + ([base] if base else []))
                for alias in node.names:
                    candidate = f"{base}.{alias.name}" if base else alias.name
                    target = resolve(candidate) or resolve(base)
                    if target:
                        add_edge(target, node, "from_import", [alias.asname or alias.name], self.function_depth > 0, self.type_only)
                    else:
                        facts[source]["external_imports"].append(base)
                        if base == "ystwin" or base.startswith("ystwin."):
                            dynamic_gaps.append({"path": source, "line": node.lineno, "kind": "unresolved_local_import", "module": candidate})

            def visit_Call(self, node):
                function = _dotted(node.func)
                if function in {"importlib.import_module", "import_module", "__import__"}:
                    value = static_path(node.args[0], assignments, source) if node.args else None
                    target = resolve(value) if value and "*" not in value else None
                    if target:
                        add_edge(target, node, "literal_dynamic_import", [], self.function_depth > 0, self.type_only)
                    elif value and "*" not in value and not value.startswith("ystwin"):
                        facts[source]["external_imports"].append(value)
                    else:
                        dynamic_gaps.append({"path": source, "line": node.lineno, "kind": "dynamic_module_import", "literal_module": value})
                if function.endswith("spec_from_file_location"):
                    value = static_path(node.args[1], assignments, source) if len(node.args) > 1 else None
                    if value in sources:
                        add_edge(value, node, "literal_file_import", [], self.function_depth > 0, self.type_only)
                    else:
                        dynamic_gaps.append({"path": source, "line": node.lineno, "kind": "dynamic_file_import", "literal_path": value})
                argument = None
                if isinstance(node.func, ast.Attribute):
                    if node.func.attr in {"to_csv", "to_json", "savefig"} and node.args:
                        argument = node.args[0]
                    elif node.func.attr in {"write_text", "write_bytes"}:
                        argument = node.func.value
                if argument is not None and source.startswith("scripts/") and source not in MECHANICAL_REFERENCES:
                    pattern = static_path(argument, assignments, source)
                    if pattern and pattern.startswith(("outputs/", "figures/", "data/")):
                        basename = PurePosixPath(pattern).name
                        if basename not in {"*", "*.csv", "*.json", "*.png", "*.svg"}:
                            for path in sorted(all_paths):
                                if fnmatch.fnmatchcase(path, pattern):
                                    writers.append({"source": source, "target": path, "line": node.lineno, "kind": "static_write_destination", "pattern": pattern})
                self.generic_visit(node)

        Visitor().visit(tree)
        if source not in MECHANICAL_REFERENCES:
            for node in ast.walk(tree):
                if isinstance(node, (ast.Constant, ast.BinOp, ast.Call)):
                    value = static_path(node, assignments, source)
                    if value in all_paths and value != source:
                        references.append({"source": source, "target": value, "line": node.lineno, "kind": "literal_path_reference_not_execution_proof"})
        facts[source]["external_imports"] = sorted(set(facts[source]["external_imports"]))
    grouped_edges = {}
    for edge in edges:
        key = tuple(edge[field] for field in ("source", "target", "kind", "line", "deferred", "type_checking_only"))
        if edge["kind"] == "package_initialization":
            key = tuple(edge[field] for field in ("source", "target", "kind", "deferred", "type_checking_only"))
        if key not in grouped_edges:
            grouped_edges[key] = deepcopy(edge)
        else:
            previous = grouped_edges[key]
            previous["bound_names"] = sorted(set(previous["bound_names"]) | set(edge["bound_names"]))
            previous["name_use_lines"] = sorted(set(previous["name_use_lines"]) | set(edge["name_use_lines"]))
            previous["line"] = min(previous["line"], edge["line"])
    unique_references = {canonical(edge): edge for edge in references}
    unique_writers = {canonical(edge): edge for edge in writers}
    return {
        "python_files": facts,
        "imports": sorted(grouped_edges.values(), key=lambda edge: (edge["source"], edge["line"], edge["target"], edge["kind"])),
        "path_references": sorted(unique_references.values(), key=lambda edge: (edge["source"], edge["line"], edge["target"])),
        "static_writers": sorted(unique_writers.values(), key=lambda edge: (edge["target"], edge["source"], edge["line"])),
        "dynamic_resolution_gaps": dynamic_gaps,
    }


def structured_imports(root: Path, contracts: list[dict], sources: dict[str, str], paths: set[str], source_hashes: dict) -> list[dict]:
    edges = []
    modules = {module_name(name): name for name in sources}
    for contract in contracts:
        table, consumer = relative_path(contract["table"]), relative_path(contract["consumer"])
        if table not in paths or consumer not in sources:
            raise InventoryError(f"Structured import contract lacks its table or consumer: {table}")
        tree = ast.parse(sources[consumer])
        functions = [node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == contract["function"]]
        calls = []
        for function in functions:
            for node in ast.walk(function):
                if not isinstance(node, ast.Call) or _dotted(node.func) != "importlib.import_module" or not node.args:
                    continue
                value = node.args[0]
                if (isinstance(value, ast.JoinedStr) and len(value.values) == 2
                        and isinstance(value.values[0], ast.Constant)
                        and value.values[0].value == contract["module_prefix"]
                        and isinstance(value.values[1], ast.FormattedValue)
                        and isinstance(value.values[1].value, ast.Attribute)
                        and value.values[1].value.attr == contract["column"]):
                    calls.append(node)
        if len(calls) != 1:
            raise InventoryError(f"Structured import expression no longer matches reviewed consumer: {consumer}")
        data = read_public(root, table)
        source_hashes[table] = hashlib.sha256(data).hexdigest()
        lines = [(index, line) for index, line in enumerate(data.decode("utf-8-sig").splitlines(), 1) if line.strip() and not line.startswith("#")]
        reader = csv.DictReader((line for _, line in lines), delimiter="\t")
        if contract["column"] not in (reader.fieldnames or []):
            raise InventoryError(f"Structured import table lacks column {contract['column']}: {table}")
        targets = defaultdict(list)
        for row in reader:
            name = row[contract["column"]]
            if name in {"", "-"}:
                continue
            module = contract["module_prefix"] + name
            target = modules.get(module)
            if target is None:
                raise InventoryError(f"Structured import table names an uninventoried module: {module}")
            targets[target].append(lines[reader.line_num - 1][0])
        for target, rows in sorted(targets.items()):
            if target != consumer:
                edges.append({
                    "source": consumer, "target": target, "kind": "data_declared_validation_import",
                    "line": calls[0].lineno, "bound_names": [], "name_use_lines": [],
                    "deferred": True, "type_checking_only": False,
                    "registry_path": table, "registry_rows": rows,
                    "use_scope": "Imports each declared module to validate registry coverage, not proof of biological execution",
                })
    return edges


def document_references(documents: dict[str, str], paths: set[str]) -> list[dict]:
    aliases = defaultdict(set)
    for name in paths:
        aliases[name].add(name)
        if name.startswith("src/ystwin/"):
            aliases[name.removeprefix("src/ystwin/")].add(name)
            dotted = module_name(name)
            if dotted:
                aliases[dotted].add(name)
    references = []
    for source, text in sorted(documents.items()):
        for number, line in enumerate(text.splitlines(), 1):
            values = [match.group(1).split("::", 1)[0] for match in re.finditer(r"`([^`\n]+)`", line)]
            values += [match.group(1).split("#", 1)[0] for match in re.finditer(r"\]\(([^)\s]+)\)", line)]
            for value in values:
                candidates = set(aliases.get(value, ()))
                relative = str(PurePosixPath(source).parent / value)
                if relative in paths:
                    candidates.add(relative)
                for target in sorted(candidates):
                    if target != source:
                        references.append({"source": source, "target": target, "line": number, "kind": "document_reference_not_execution_proof"})
    return references


def reachable(roots: list[str], imports: list[dict], include_validation_imports=True) -> set[str]:
    adjacency = defaultdict(set)
    for edge in imports:
        validation_only = edge["kind"] == "data_declared_validation_import"
        if not edge["type_checking_only"] and (include_validation_imports or not validation_only):
            adjacency[edge["source"]].add(edge["target"])
    result, pending = set(), list(roots)
    while pending:
        name = pending.pop()
        if name not in result:
            result.add(name)
            pending.extend(adjacency[name] - result)
    return result


def strongly_connected(imports: list[dict]) -> list[list[str]]:
    graph = defaultdict(set)
    for edge in imports:
        if edge["source"].startswith("src/ystwin/") and edge["target"].startswith("src/ystwin/") and not edge["type_checking_only"]:
            graph[edge["source"]].add(edge["target"])
            graph[edge["target"]]
    indices, lowlinks, stack, active, groups = {}, {}, [], set(), []

    def visit(name):
        indices[name] = lowlinks[name] = len(indices)
        stack.append(name)
        active.add(name)
        for target in sorted(graph[name]):
            if target not in indices:
                visit(target)
                lowlinks[name] = min(lowlinks[name], lowlinks[target])
            elif target in active:
                lowlinks[name] = min(lowlinks[name], indices[target])
        if indices[name] == lowlinks[name]:
            component = []
            while True:
                target = stack.pop()
                active.remove(target)
                component.append(target)
                if target == name:
                    break
            if len(component) > 1:
                groups.append(sorted(component))

    for name in sorted(graph):
        if name not in indices:
            visit(name)
    return sorted(groups)


def validate_workflow(workflow: dict, paths: set[str]) -> tuple[list[str], dict, dict]:
    if workflow.get("schema_version") != SCHEMA_VERSION or workflow.get("scope_complete") is not False:
        raise InventoryError("Workflow must retain the current schema and scope_complete=false")
    nodes = workflow.get("nodes", [])
    ids = [node["id"] for node in nodes]
    if len(ids) != len(set(ids)) or not ids:
        raise InventoryError("Workflow node IDs must be unique and nonempty")
    owners, membership = {}, defaultdict(set)
    all_ids = set(ids)
    for node in nodes:
        for field in ("deliverable", "owner", "status", "completion_conditions", "tests"):
            if not node.get(field):
                raise InventoryError(f"Workflow {node['id']} is missing {field}")
        if node["status"] == "completed" and not node.get("acceptance_evidence"):
            raise InventoryError(f"Workflow {node['id']} cannot be complete without acceptance evidence")
        if set(node.get("depends_on", [])) - all_ids:
            raise InventoryError(f"Unknown dependency in {node['id']}")
        owners[node["id"]] = node["owner"]
        for name in node.get("existing_paths", []) + node["tests"]:
            relative_path(name)
            if name not in paths:
                raise InventoryError(f"Workflow references missing path: {name}")
        for name in node.get("existing_paths", []) + node["tests"]:
            membership[name].add(node["id"])
        for module in node.get("existing_modules", []):
            name = "src/" + module.replace(".", "/") + ".py"
            if name not in paths:
                raise InventoryError(f"Workflow references missing module: {module}")
            membership[name].add(node["id"])
    additions = set()
    for addition in workflow.get("reviewed_session_additions", []):
        name = relative_path(addition["path"])
        if name in additions:
            raise InventoryError(f"Duplicate reviewed session addition: {name}")
        additions.add(name)
        if name not in paths or addition["work_item"] not in all_ids:
            raise InventoryError(f"Unknown reviewed session addition or work owner: {name}")
        if "purpose" in addition and (not isinstance(addition["purpose"], str) or not addition["purpose"].strip()):
            raise InventoryError(f"Reviewed session addition lacks its use contract: {name}")
        membership[name].add(addition["work_item"])
    try:
        order = list(TopologicalSorter({node["id"]: node.get("depends_on", []) for node in nodes}).static_order())
    except CycleError as exc:
        raise InventoryError(f"Implementation workflow has a dependency cycle: {exc}") from exc
    return order, owners, membership


def cleanup_command(name: str) -> list[str]:
    contract = CLEANUP_CONTRACTS[name]
    if contract["action"] == "relocate_immutable_source_evidence":
        return ["git", "mv", "--", name, contract["replacement_paths"][0]]
    return ["git", "rm", "--", name] if contract["original"]["tracked"] else ["rm", "--", name]


def verify_cleanup_bytes(data: bytes, identity: dict, name: str) -> None:
    if len(data) != identity["size_bytes"] or hashlib.sha256(data).hexdigest() != identity["sha256"]:
        raise InventoryError(f"Approved cleanup original/source identity mismatch: {name}")


def validate_task_replacement(root: Path, receipt: dict) -> None:
    data = read_public(root, TASK_REPORT)
    if receipt.get("replacement_sha256") != hashlib.sha256(data).hexdigest():
        raise InventoryError("Task replacement is not bound to its execution receipt")
    report = parse_json(data)
    tasks = ["biomass", "ethanol", "CO2"]
    models = report.get("models", {})
    if (report.get("tasks") != tasks or report.get("status") != "completed"
            or report.get("biological_validation") is not False or set(models) != {"plain", "ec"}):
        raise InventoryError("Task replacement must retain the completed matched-host full-task refusal")
    for model in models.values():
        if (model.get("tasks") != tasks or model.get("status") != "infeasible"
                or model.get("order") != [] or not model.get("refusal")
                or model.get("constrained_growth") is not None or model.get("biological_validation") is not False):
            raise InventoryError("Task replacement lost a measured task or its infeasibility finding")
    claims = parse_json(read_public(root, "data/current_claims.json"))
    evidence = claims.get("profiles", {}).get("cosmic", {}).get("evidence", [])
    if TASK_REPORT not in {item.get("path") for item in evidence} or TASK_TABLE in {item.get("path") for item in evidence}:
        raise InventoryError("Current COSMIC claims have not adopted the replacement report")
    registry = parse_json(read_public(root, "data/artifact_registry.json"))
    artifacts = [item for run in registry.get("runs", []) for item in run.get("artifacts", [])]
    if TASK_REPORT not in {item.get("path") for item in artifacts}:
        raise InventoryError("Artifact provenance has not registered the replacement report")


def validate_cleanup(root: Path, document: dict, tracked: dict, baseline: dict) -> dict:
    if document.get("approval_required") is not True or document.get("execution_allowed") is not False:
        raise InventoryError("Deletion manifest must keep automatic execution disabled and require approval")
    result = {"approved": {}, "completed": {}, "relocations": {}, "replacements": {}}
    requests = [r for r in [document.get("approval_request")] if r is not None]
    requests += document.get("additional_approval_requests", [])
    if not requests:
        return result
    for request in requests:
        _validate_cleanup_request(root, request, tracked, baseline, result)
    if set(result["approved"]) != set(CLEANUP_CONTRACTS):
        raise InventoryError(
            f"Cleanup requests must retain all {len(CLEANUP_CONTRACTS)} exact approved actions")
    executed = [name for request in requests for name in request.get("executed_paths", [])]
    if len(executed) != len(set(executed)) or set(executed) != set(result["completed"]):
        raise InventoryError("Executed cleanup paths must match distinct approved execution receipts exactly")
    return result


def _validate_cleanup_request(root: Path, request: dict, tracked: dict, baseline: dict, result: dict) -> None:
    """Validate one approval batch. Each contract names its own batch, so a later approval
    cannot be backdated into an earlier batch's sealed record."""
    approval = CLEANUP_APPROVALS.get(request.get("id"))
    if approval is None or request.get("status") != "approved" or canonical(request.get("approval")) != canonical(approval):
        raise InventoryError("Unknown or unapproved exact cleanup request")
    batch_id = request["id"]
    batch = {name for name, c in CLEANUP_CONTRACTS.items() if c["batch"] == batch_id}
    approved_here = set()
    for action in request.get("actions", []):
        name = relative_path(action["path"])
        contract = CLEANUP_CONTRACTS.get(name)
        if (contract is None or name in result["approved"] or action.get("action") != contract["action"]
                or action.get("replacement_paths") != contract["replacement_paths"]):
            raise InventoryError(f"Cleanup approval action/replacement does not match its exact reviewed scope: {name}")
        if name == RELEASE_SOURCE and (action.get("replacement") != RELEASE_ASSET
                                      or action.get("sha256") != contract["original"]["sha256"]):
            raise InventoryError("Immutable source relocation differs from its approved path/hash")
        result["approved"][name] = action
        approved_here.add(name)
    if approved_here != batch:
        raise InventoryError(
            f"Cleanup request {batch_id} must retain all {len(batch)} of its exact approved actions")
    for receipt in request.get("execution_receipts", []):
        name = relative_path(receipt["path"])
        contract = CLEANUP_CONTRACTS.get(name)
        if contract is None or name in result["completed"]:
            raise InventoryError(f"Unknown or duplicate cleanup execution receipt: {name}")
        if (receipt.get("approval_request_id") != batch_id or contract["batch"] != batch_id
                or receipt.get("action") != contract["action"]
                or canonical(receipt.get("original")) != canonical(contract["original"])
                or receipt.get("replacement_paths") != contract["replacement_paths"]
                or receipt.get("command") != cleanup_command(name)
                or type(receipt.get("returncode")) is not int or receipt["returncode"] != 0
                or not receipt.get("operator") or receipt.get("date") != approval["date"]):
            raise InventoryError(f"Cleanup execution receipt does not match the approved exact action: {name}")
        path = root / name
        if name in tracked or path.exists() or path.is_symlink() or path.resolve() != path:
            raise InventoryError(f"Completed cleanup path is still tracked, present or aliased: {name}")
        identity = contract["original"]
        if identity["tracked"]:
            original = baseline.get(name, {})
            if original.get("kind") != "blob" or original.get("object_id") != identity["git_object_id"]:
                raise InventoryError(f"Cleanup original Git identity is not the reviewed baseline: {name}")
            data = git(root, "show", f"{identity['git_commit']}:{name}")
            verify_cleanup_bytes(data, identity, name)
        for replacement in contract["replacement_paths"]:
            safe_file(root, replacement)
        if contract["action"] == "relocate_immutable_source_evidence":
            replacement = contract["replacement_paths"][0]
            if (tracked.get(replacement, {}).get("object_id") != identity["git_object_id"]
                    or tracked.get(replacement, {}).get("mode") != baseline[name].get("mode")):
                raise InventoryError(f"Relocated source is not tracked with its original Git blob/mode: {replacement}")
            verify_cleanup_bytes(read_public(root, replacement), identity, replacement)
            result["relocations"][name] = replacement
        elif contract["action"] == "replace_superseded_result":
            validate_task_replacement(root, receipt)
            result["replacements"][name] = contract["replacement_paths"][0]
        result["completed"][name] = deepcopy(receipt)


def frozen_evidence(root: Path, paths: set[str], baseline: dict, relocations=None) -> dict:
    relocations = relocations or {}
    if any(name != RELEASE_SOURCE or asset != RELEASE_ASSET for name, asset in relocations.items()):
        raise InventoryError("Unreviewed historical source relocation")
    required = {name for name in baseline if name.startswith(("archive/", "data/frozen_evidence/native_v1/")) or (name.endswith(".py") and context_kind(name) == "historical_evidence")}
    missing = {name for name in required if relocations.get(name, name) not in paths}
    if missing:
        raise InventoryError(f"Protected historical paths are no longer inventoried: {sorted(missing)}")
    records, references, checks = {}, defaultdict(list), []
    manifest_path = FROZEN_MANIFEST if FROZEN_MANIFEST in paths else None
    if manifest_path:
        data = read_public(root, manifest_path)
        if hashlib.sha256(data).hexdigest() != FROZEN_MANIFEST_SHA256:
            raise InventoryError("Public frozen manifest does not match its independently pinned SHA-256")
        manifest = parse_json(data)
        directory = str(PurePosixPath(manifest_path).parent)
        checks.append({"path": manifest_path, "sha256": FROZEN_MANIFEST_SHA256, "basis": "independently_pinned_public_manifest", "verified": True})
        for record in manifest["records"]:
            public_path = relative_path(record["public_path"])
            name = relative_path(directory + "/" + public_path)
            if name not in paths:
                raise InventoryError(f"Frozen public artifact is not inventoried: {name}")
            if name in records:
                raise InventoryError(f"Duplicate frozen public artifact: {name}")
            artifact = read_public(root, name)
            actual = hashlib.sha256(artifact).hexdigest()
            if actual != record["public_sha256"] or len(artifact) != record["public_size_bytes"]:
                raise InventoryError(f"Frozen public artifact changed: {name}")
            records[name] = {
                "origin_path": record["origin_path"],
                "original_sha256": record["original_sha256"],
                "public_sha256": record["public_sha256"],
                "format": record["format"],
            }
            references[record["origin_path"]].append(name)
            checks.append({"path": name, "sha256": actual, "basis": "public_manifest_record", "verified": True})
        for reference in manifest["external_references"]:
            references[reference["origin_path"]].append("lineage_only_boundary_not_read")
    for name in sorted(paths | set(relocations)):
        historical = name.startswith("archive/") or (name.endswith(".py") and context_kind(name) == "historical_evidence")
        if historical and name in baseline and baseline[name]["kind"] == "blob":
            storage = relocations.get(name, name)
            data = read_public(root, storage)
            algorithm = hashlib.sha1 if len(baseline[name]["object_id"]) == 40 else hashlib.sha256
            blob = algorithm(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            if blob != baseline[name]["object_id"]:
                raise InventoryError(f"Historical archive/executor changed from baseline: {name}")
            check = {"path": storage, "sha256": hashlib.sha256(data).hexdigest(), "basis": "baseline_historical_blob", "verified": True}
            if storage != name:
                original = root / name
                if name in paths or original.exists() or original.is_symlink() or original.resolve() != original:
                    raise InventoryError(f"Ambiguous or aliased historical source relocation: {name}")
                verify_cleanup_bytes(data, CLEANUP_CONTRACTS[name]["original"], storage)
                check.update(logical_source_path=name, approval_request_id=CLEANUP_APPROVAL_ID,
                             original_git_commit=CLEANUP_BASELINE, original_git_object_id=blob,
                             basis="approved_byte_preserving_relocation_of_baseline_historical_blob")
            checks.append(check)
    return {"manifest": manifest_path, "records": records, "references": dict(references), "integrity_checks": checks}


def check_frozen_claims(root: Path, assertions: list[dict], frozen: dict) -> list[dict]:
    checks, payloads = [], {}
    for assertion in assertions:
        name = relative_path(assertion["path"])
        if name not in frozen["records"]:
            raise InventoryError(f"Negative-claim assertion is not bound to pinned public evidence: {name}")
        if name not in payloads:
            data = read_public(root, name)
            if hashlib.sha256(data).hexdigest() != frozen["records"][name]["public_sha256"]:
                raise InventoryError(f"Pinned claim artifact changed during the scan: {name}")
            payloads[name] = parse_json(data)
        pointer = assertion["pointer"]
        if not isinstance(pointer, str) or not pointer.startswith("/"):
            raise InventoryError("Frozen claim assertion requires a JSON pointer")
        value = payloads[name]
        try:
            for component in pointer[1:].split("/"):
                key = component.replace("~1", "/").replace("~0", "~")
                if isinstance(value, list):
                    if not key.isdigit():
                        raise KeyError(key)
                    value = value[int(key)]
                else:
                    value = value[key]
        except (KeyError, IndexError, TypeError) as exc:
            raise InventoryError(f"Frozen claim pointer is absent: {name}{pointer}") from exc
        if canonical(value) != canonical(assertion["equals"]):
            raise InventoryError(f"Retained negative claim contradicts pinned evidence: {name}{pointer}")
        checks.append({**assertion, "verified": True})
    return checks


def expand_negative_registry(registry: list[dict], paths: set[str], baseline_paths: set[str], frozen: dict, root: Path, relocations=None, replacements=None) -> tuple[list[dict], dict]:
    relocations = relocations or {}
    replacements = replacements or {}
    storage_paths = {**relocations, **replacements}
    expanded, memberships, ids = [], defaultdict(set), set()
    for original in registry:
        item = deepcopy(original)
        identity = item.get("id")
        if not identity or identity in ids:
            raise InventoryError("Retained-negative IDs must be unique and nonempty")
        ids.add(identity)
        for field in ("claim", "scope", "supporting_tests", "claim_documents", "preservation"):
            if not item.get(field):
                raise InventoryError(f"Negative {identity} lacks {field}")
        explicit = set(item.get("artifact_paths", []))
        historical = {relative_path(name) for name in item.get("historical_artifact_paths", [])}
        for name in historical:
            if name not in replacements or replacements[name] not in explicit:
                raise InventoryError(f"Negative {identity} historical artifact lacks an approved completed replacement retained in artifact_paths: {name}")
        for name in explicit | set(item.get("supporting_code", [])) | set(item["supporting_tests"]) | set(item["claim_documents"]):
            relative_path(name)
            if storage_paths.get(name, name) not in paths:
                raise InventoryError(f"Negative {identity} references missing path: {name}")
        for prefix in item.get("artifact_prefixes", []):
            relative_path(prefix.rstrip("/"))
            matched = {name for name in baseline_paths if name.startswith(prefix)}
            if not matched:
                raise InventoryError(f"Negative {identity} has an empty baseline selector: {prefix}")
            explicit.update(matched)
        for name, record in frozen["records"].items():
            if record["origin_path"] in explicit:
                explicit.add(name)
        item["retained_logical_artifact_paths"] = sorted(explicit | historical)
        item["source_asset_relocations"] = [
            {"logical_source_path": name, "storage_path": relocations[name],
             "approval_request_id": CLEANUP_APPROVAL_ID, **CLEANUP_CONTRACTS[name]["original"]}
            for name in sorted(explicit & relocations.keys())
        ]
        item["approved_result_replacements"] = [
            {"original_path": name, "replacement_path": replacements[name],
             "approval_request_id": CLEANUP_APPROVAL_ID, **CLEANUP_CONTRACTS[name]["original"],
             "execution_receipt_manifest": "data/deletion_manifest.json",
             "identity_scope": "original Git evidence plus the approved adopted result, not restoration of obsolete ranking"}
            for name in sorted((explicit | historical) & replacements.keys())
        ]
        explicit = {storage_paths.get(name, name) for name in explicit}
        if explicit - paths:
            raise InventoryError(f"Negative {identity} lost retained baseline paths: {sorted(explicit - paths)}")
        item["retained_artifact_paths"] = sorted(explicit)
        item["frozen_claim_checks"] = check_frozen_claims(root, item.get("frozen_assertions", []), frozen)
        if any(assertion["path"] not in explicit for assertion in item.get("frozen_assertions", [])):
            raise InventoryError(f"Negative {identity} asserts an artifact it does not retain")
        item["approval_required_for_any_retirement"] = True
        item["tests_are_scoped_software_checks_not_biological_authorization"] = True
        for name in explicit:
            memberships[name].add(identity)
        expanded.append(item)
    if TASK_TABLE in replacements:
        n06 = next((item for item in expanded if item["id"] == "N06"), {})
        if (TASK_REPORT not in n06.get("artifact_paths", [])
                or TASK_TABLE in n06.get("artifact_paths", [])
                or TASK_TABLE not in n06.get("historical_artifact_paths", [])):
            raise InventoryError("Approved task-result replacement must preserve N06 with the adopted report and original historical identity")
    return expanded, memberships


def semantic_graph(imports: list[dict]) -> list[dict]:
    rows = {}
    for edge in imports:
        row = {key: edge[key] for key in ("source", "target", "kind", "bound_names", "deferred", "type_checking_only")}
        rows[canonical(row)] = row
    return sorted(rows.values(), key=canonical)


def file_kind(name: str) -> str:
    if name.endswith(".py"):
        return "test" if context_kind(name) == "test" or PurePosixPath(name).name.startswith("test_") else "python_source"
    if name.endswith((".md", ".Rd")):
        return "document"
    if name.startswith("src/ystwin.egg-info/"):
        return "generated_package_metadata"
    if name.startswith(("outputs/", "figures/", "archive/outputs/")):
        return "committed_result"
    if name.endswith((".m", ".mltbx")):
        return "external_language_source"
    if name.startswith("data/"):
        return "data_or_evidence"
    return "configuration_or_other"


def classify(name: str, policy: dict, membership: dict, active: set[str], validation: set[str], negative: dict, frozen: dict) -> tuple[str, str, set[str]]:
    work = set(membership.get(name, ()))
    for explicit in policy["explicit_path_roles"]:
        if name == explicit["path"]:
            return explicit["role"], explicit["purpose"], work | {explicit["work_item"]}
    if name.startswith("src/ystwin.egg-info/"):
        return "superseded", "Generated setuptools metadata; pyproject/package discovery is authoritative. Retain until approval and clean-build equivalence review.", {"W01"}
    if name.startswith(("archive/code/", "archive/tests/", "archive/docs/", "docs/superseded/")):
        return "superseded", "Explicit historical retirement/correction record, preserved rather than imported or treated as a current assertion.", {"W13", "W15"}
    if name in negative:
        return "recorded_negative_result", "Retained evidence/support for " + ", ".join(sorted(negative[name])) + "; includes failed applicability and refusal, not a successful biological claim.", work | {"W13"}
    if name.startswith("src/ystwin/") and name.endswith("/__init__.py"):
        if name in active:
            return "active_runtime", "Package initialization/re-export on a supported runtime import path.", work or {"W09"}
        return "supported_public_library", "Package namespace required for retained public submodules; namespace presence is not a completed scientific layer.", work or {"W14"}
    if name in active:
        return "active_runtime", "Static dependency of an inspected supported runtime root; optional/deferred reachability does not assert that a branch executes or sets a prediction.", work
    for public in policy["supported_public_apis"]:
        if public["path"] == name:
            return "supported_public_library", public["purpose"], work | {public["work_item"]}
    if name.startswith("tests/"):
        return "validation_reproduction", "Live test or fixture infrastructure. Its imports are validation evidence, not reasons to retain otherwise unsupported runtime modules.", work or {"W15"}
    if name in validation and name.startswith(("src/ystwin/", "scripts/")):
        return "validation_reproduction", "Static dependency of a reviewed calibration, audit, simulation, scoring or reproduction runner, separate from the runtime roots.", work
    if name.startswith("scripts/") and work:
        return "validation_reproduction", "Reviewed operational/provisioning/reproduction script with an explicit work owner; not automatically part of the product runtime.", work
    if name.startswith("src/ystwin/"):
        return "unowned_requires_decision", "No verified non-test root or declared public-extension use contract establishes this module's current role. Retain pending owner review; this is not an orphan/deletion finding.", work or {"W15"}
    for domain in policy["reviewed_data_domains"]:
        if name.startswith(domain["prefix"]):
            role = "active_runtime" if domain["prefix"] == "data/pathways/" else "validation_reproduction"
            return role, domain["purpose"] + "; exact producer/input provenance is a separate artifact-registry responsibility.", work | {domain["work_item"]}
    if name.startswith(("outputs/", "figures/")):
        return "validation_reproduction", "Committed result retained for comparison/reproduction. Static writer evidence is reported separately; this role does not assert a current rerun or complete provenance.", work or {"W13", "W15"}
    if name.endswith(".md"):
        return "validation_reproduction", "Operational, scientific or historical documentation. Located references are recorded; prose correctness/completion is not inferred from file existence.", work or {"W15"}
    if name in OWNED_PATHS:
        return "validation_reproduction", "Inventory policy, work DAG or approval-only retirement review owned by W00.", {"W00"}
    if name.startswith("data/") and work:
        return "validation_reproduction", "Explicitly reviewed session registry/evidence metadata with an accountable owner. Its content and any scientific acceptance remain that owner's responsibility.", work
    if name.startswith(".github/") or name in {"pyproject.toml", "requirements-quality.txt", ".gitignore", ".ystwin.env.example"}:
        return "validation_reproduction", "Declared installation, public provisioning or CI configuration; no secret configuration values are inspected.", work or {"W01"}
    return "unowned_requires_decision", "No reviewed role for this tracked path. Retain and assign an explicit purpose/owner before any migration decision.", work or {"W15"}


def validate_deletions(document: dict, paths: set[str], protected: set[str], graph: dict, frozen: dict, negatives: dict, cleanup: dict) -> dict:
    result = deepcopy(document)
    seen = set()
    references = graph["path_references"] + graph.get("document_references", [])
    for proposal in result.get("proposals", []):
        name = relative_path(proposal["path"])
        completed = name in cleanup["completed"]
        approved = name in cleanup["approved"]
        if any(char in name for char in "*?[") or (name not in paths and not completed) or name in seen:
            raise InventoryError(f"Deletion proposal must name one unique inventoried or approved-completed exact path: {name}")
        seen.add(name)
        if (proposal.get("approval_required") is not True or proposal.get("approved") is not approved
                or proposal.get("action_executed", False) is not completed):
            raise InventoryError(f"Deletion proposal cannot forge approval or execution: {name}")
        if not proposal.get("reason") or not proposal.get("replacement_contract") or not proposal.get("approval_conditions"):
            raise InventoryError(f"Deletion proposal needs reason, replacement and approval conditions: {name}")
        if approved:
            contract = CLEANUP_CONTRACTS[name]
            if (proposal.get("approval_request_id") != contract["batch"]
                    or proposal.get("proposed_action") != contract["action"]
                    or proposal.get("replacement_paths") != contract["replacement_paths"]):
                raise InventoryError(f"Deletion proposal differs from its exact approved action: {name}")
        else:
            for replacement in proposal.get("replacement_paths", []):
                if relative_path(replacement) not in paths:
                    raise InventoryError(f"Unknown deletion replacement: {replacement}")
        if name in protected and not (approved and name in {RELEASE_SOURCE, TASK_TABLE}):
            raise InventoryError(f"Deletion proposal targets protected historical/negative evidence: {name}")
        storage = cleanup["relocations"].get(name, cleanup["replacements"].get(name, name))
        callers = [edge for edge in graph["imports"] if edge["target"] == name]
        proposal["callers"] = {
            "static_imports": callers,
            "literal_or_document_references": [edge for edge in references if edge["target"] == name and edge["source"] not in MECHANICAL_REFERENCES],
            "interpretation": "References include logical historical identities, not only runtime use. Absence of static callers never supplies approval; exact user review and replacement contracts remain required.",
        }
        proposal["frozen_impact"] = {
            "public_bundle_or_lineage_references": frozen["references"].get(name, []),
            "retained_negative_ids": sorted(negatives.get(storage, ())),
            "is_protected_evidence": name in protected or storage in protected,
            "assessment": "Only the independently recognized user-approved exact action can change storage. Historical logical identities and negative findings remain retained; no current-source freshness exemption is created.",
        }
        proposal["action_executed"] = completed
    if cleanup["approved"] and seen != set(cleanup["approved"]):
        raise InventoryError("Every approved exact cleanup action requires its own matching proposal record")
    if "approval_request" in result:
        result["approval_request"]["pending_paths"] = sorted(set(cleanup["approved"]) - set(cleanup["completed"]))
    result["execution_verification"] = [
        {"path": name, "status": "verified_completed", "original_sha256": receipt["original"]["sha256"],
         "original_absent_and_untracked": True, "replacement_files_present": True}
        for name, receipt in sorted(cleanup["completed"].items())
    ]
    result["protected_exact_paths"] = sorted(protected)
    return result


def build_inventory(root: Path, seed: dict, workflow: dict, deletions: dict, include_untracked=(), snapshot_hashes=False) -> tuple[dict, dict, dict]:
    root = root.resolve()
    policy = seed["classification_policy"]
    if tuple(policy.get("roles", ())) != ROLES:
        raise InventoryError("Inventory role vocabulary changed or is incomplete")
    tracked, baseline, scan = git_catalog(root, workflow["baseline"]["commit"])
    if deletions.get("approval_request") and workflow["baseline"]["commit"] != CLEANUP_BASELINE:
        raise InventoryError("Reviewed cleanup must retain the original baseline commit")
    cleanup = validate_cleanup(root, deletions, tracked, baseline)
    completed = set(cleanup["completed"])
    relocations = cleanup["relocations"]
    paths = (set(tracked) | set(OWNED_PATHS) | {"S.csv", "jws.html"}) - completed
    reviewed_additions = {relative_path(item["path"]) for item in workflow.get("reviewed_session_additions", [])}
    extensions = sorted({relative_path(name) for name in include_untracked})
    for name in set(extensions) | reviewed_additions:
        if not python_allowed(name) and not name.startswith((".github/", "data/", "outputs/")) and name != "requirements-quality.txt":
            raise InventoryError(f"Untracked scan inclusion must be explicitly public code/CI or artifact metadata: {name}")
        safe_file(root, name)
        paths.add(name)
    special = [path_status(root, name, set(tracked)) for name in ("S.csv", "jws.html")]
    missing = sorted(name for name in tracked if not (root / name).exists() and not (root / name).is_symlink())
    if missing:
        raise InventoryError(f"Tracked paths missing from worktree without completed exact cleanup: {missing}")
    if set(baseline) - set(tracked) - completed:
        raise InventoryError(f"Baseline paths lost without approved completed exact actions: {sorted(set(baseline) - set(tracked) - completed)}")
    order, owners, membership = validate_workflow(workflow, paths)
    decision_ids = set()
    for decision in seed.get("decisions", []):
        if not decision.get("id") or decision["id"] in decision_ids or not decision.get("owner") or not decision.get("question"):
            raise InventoryError("Decisions require unique IDs, accountable owners and explicit questions")
        decision_ids.add(decision["id"])
        if decision.get("status") not in {"open", "resolved", "deferred"} or set(decision.get("blocks", [])) - set(owners):
            raise InventoryError(f"Invalid decision state or work dependency: {decision['id']}")
        for name in decision.get("paths", []):
            if relative_path(name) not in paths and name not in completed:
                raise InventoryError(f"Decision {decision['id']} refers to an uninventoried path: {name}")
            if name in completed and (decision["status"] != "resolved" or decision.get("approval_request_id") != CLEANUP_APPROVAL_ID):
                raise InventoryError(f"Completed cleanup decision lacks its resolved approval: {name}")
    sources, documents, source_hashes = {}, {}, {}
    for name in sorted(paths):
        if name in set(policy.get("do_not_read", ())) | STATUS_ONLY:
            continue
        if python_allowed(name) or name.endswith(".md"):
            data = read_public(root, name)
            source_hashes[name] = hashlib.sha256(data).hexdigest()
            try:
                text = data.decode("utf-8-sig")
            except UnicodeDecodeError as exc:
                raise InventoryError(f"Public source is not UTF-8: {name}") from exc
            if python_allowed(name):
                sources[name] = text
            else:
                documents[name] = text
    graph = analyze_python(sources, paths | completed)
    declared_imports = structured_imports(root, policy.get("structured_imports", []), sources, paths, source_hashes)
    graph["imports"].extend(declared_imports)
    resolved_calls = {(edge["source"], edge["line"]) for edge in declared_imports}
    graph["dynamic_resolution_gaps"] = [gap for gap in graph["dynamic_resolution_gaps"] if (gap["path"], gap["line"]) not in resolved_calls]
    for gap in graph["dynamic_resolution_gaps"]:
        gap["context"] = context_kind(gap["path"])
        gap["scope"] = "historical non-collected import, not a live missing dependency" if gap["context"] == "historical_evidence" else "static-resolution limitation requiring owner review, not proof of a broken runtime"
    graph["document_references"] = document_references(documents, paths | completed)
    graph["source_strongly_connected_components"] = strongly_connected(graph["imports"])
    frozen = frozen_evidence(root, paths, baseline, relocations)
    registry, negative = expand_negative_registry(seed["retained_negative_registry"], paths, set(baseline), frozen, root,
                                                 relocations, cleanup["replacements"])
    for record in registry:
        subjects = {"W13"}
        for name in record["supporting_code"]:
            subjects.update(membership.get(name, ()))
        record["associated_work_items"] = sorted(subjects)
        for name in record["retained_artifact_paths"]:
            membership[name].update(subjects)
    for edge in graph["static_writers"]:
        membership[edge["target"]].update(membership.get(edge["source"], ()))
    for edge in graph["imports"]:
        if context_kind(edge["source"]) == "test":
            membership[edge["source"]].update(membership.get(edge["target"], ()))
    roots = policy["active_roots"]
    if set(roots) - paths:
        raise InventoryError(f"Missing declared runtime roots: {sorted(set(roots) - paths)}")
    active = reachable(roots, graph["imports"], include_validation_imports=False)
    validation_roots = sorted(name for name in sources if name.startswith("scripts/") and name in membership and name not in roots)
    validation = reachable(validation_roots, graph["imports"])
    incoming = defaultdict(lambda: defaultdict(set))
    references = defaultdict(list)
    producer_map = defaultdict(list)
    for edge in graph["imports"]:
        incoming[edge["target"]][context_kind(edge["source"])].add(edge["source"])
    for edge in graph["path_references"] + graph["document_references"]:
        references[edge["target"]].append(edge)
    for edge in graph["static_writers"]:
        producer_map[edge["target"]].append(edge)
    entries = []
    protected = set(negative) | set(frozen["records"]) | set(relocations) | set(relocations.values()) | {name for name in paths if name.startswith("archive/")}
    source_identities = {
        asset: {"logical_source_path": name, "storage_path": asset, "materialized_path": name,
                "approval_request_id": CLEANUP_APPROVAL_ID, **CLEANUP_CONTRACTS[name]["original"],
                "identity_scope": "immutable original Git source, never mutable-runtime freshness"}
        for name, asset in relocations.items()
    }
    if frozen["manifest"]:
        protected.add(frozen["manifest"])
    reviewed_uses = {item["path"]: item["purpose"] for item in workflow.get("reviewed_session_additions", []) if "purpose" in item}
    for name in sorted(paths):
        role, purpose, work = classify(name, policy, membership, active, validation, negative, frozen)
        purpose = reviewed_uses.get(name, purpose)
        if role not in ROLES:
            raise InventoryError(f"Invalid role for {name}: {role}")
        if not work:
            work = {"W15"}
        if work - set(owners):
            raise InventoryError(f"Unknown work ownership for {name}: {sorted(work)}")
        entry = {
            "path": name,
            "tracked": name in tracked,
            "kind": file_kind(name),
            "role": role,
            "purpose": purpose,
            "work_items": sorted(work),
            "accountable_owners": sorted({owners[item] for item in work}),
            "baseline_git_object_id": baseline.get(name, {}).get("object_id"),
            "observed_index_object_id": tracked.get(name, {}).get("object_id"),
            "hash_scope": "Git object identity only; not a mutable-worktree freshness requirement",
            "imported_by": {kind: sorted(names) for kind, names in sorted(incoming[name].items())},
            "reference_count": len(references[name]),
            "retained_negative_ids": sorted(negative.get(name, ())),
            "retirement_protected": name in protected,
            "decision_ids": [decision["id"] for decision in seed.get("decisions", []) if name in decision.get("paths", [])],
        }
        if name in source_identities:
            entry["source_identity"] = source_identities[name]
            entry["baseline_git_object_id"] = source_identities[name]["git_object_id"]
            entry["hash_scope"] = source_identities[name]["identity_scope"]
        if name in graph["python_files"]:
            entry["python"] = graph["python_files"][name]
        if name in frozen["records"]:
            entry["frozen_record"] = frozen["records"][name]
        if file_kind(name) == "committed_result":
            entry["producer_evidence"] = producer_map[name]
            entry["reproduction_status"] = "retained_historical_or_negative_not_refreshed" if name in protected else "static_destination_only_not_a_verified_rerun" if producer_map[name] else "producer_not_established_by_this_scan_consult_artifact_registry"
        if name in STATUS_ONLY:
            entry["contents_read"] = False
        entries.append(entry)
    graph["runtime_roots"] = roots
    graph["validation_roots"] = validation_roots
    graph["limitations"] = [
        "AST imports/name loads and declared literal paths are static evidence, not traced execution or proof that optional branches run.",
        "Name uses are conservative lexical references, not a scope-resolved call graph; alias shadowing and wildcard re-exports require review.",
        "Dynamic imports with nonliteral targets and computed paths are listed as gaps, never silently used to declare a file dead.",
        "The source graph may contain deferred-import cycles. Only the implementation workflow is required to be a DAG.",
        "Static write destinations are not a complete producer registry and never certify current output freshness; artifact_registry belongs to the provenance owner.",
        "Tests are validation consumers and cannot alone justify a runtime role; public extensions require separately declared use contracts.",
    ]
    graph.pop("python_files")
    scan.update({
        "source_of_truth": "git ls-files --stage -z for tracked coverage; git ls-tree for the fixed baseline; AST of explicitly safe public worktree source for use evidence",
        "source_scan_sha256": digest(source_hashes),
        "source_scan_hash_policy": "Observed scan fingerprint only. Default --check compares import topology/roles and frozen integrity, not changing source bytes or line numbers. Use --snapshot-hashes and --check --verify-source-hashes for an exact final snapshot.",
        "source_file_count": len(source_hashes),
        "semantic_import_graph_sha256": digest(semantic_graph(graph["imports"])),
        "explicit_untracked_inclusions": extensions,
        "reviewed_session_additions": sorted(reviewed_additions),
        "owned_untracked_policy_paths": sorted(set(OWNED_PATHS) - set(tracked)),
        "untracked_discovery_policy": "Untracked files are never recursively read. The five owned files, status-only root artifacts and exact workflow.reviewed_session_additions are observed; other public code requires --include-untracked PATH. Non-source artifact inclusions are metadata-only.",
        "secrets_read": False,
        "frozen_integrity_checks": frozen["integrity_checks"],
    })
    if snapshot_hashes:
        scan["source_hashes_sha256"] = source_hashes
    result = {
        "schema_version": SCHEMA_VERSION,
        "kind": "repository_inventory",
        "classification_policy": deepcopy(policy),
        "retained_negative_registry": registry,
        "source_identity_registry": source_identities,
        "completed_cleanup_actions": [cleanup["completed"][name] for name in sorted(completed)],
        "decisions": deepcopy(seed.get("decisions", [])),
        "scan": scan,
        "special_root_artifacts": special,
        "entries": entries,
        "use_graph": graph,
        "summary": {
            "tracked_paths": len(tracked),
            "inventoried_paths": len(entries),
            "tracked_by_kind": {kind: sum(entry["tracked"] and entry["kind"] == kind for entry in entries) for kind in sorted({entry["kind"] for entry in entries})},
            "tracked_by_role": {role: sum(entry["tracked"] and entry["role"] == role for entry in entries) for role in ROLES},
            "all_by_role": {role: sum(entry["role"] == role for entry in entries) for role in ROLES},
            "live_python_modules": sum(name.startswith("src/ystwin/") and not name.endswith("/__init__.py") for name in sources),
            "live_package_initializers": sum(name.startswith("src/ystwin/") and name.endswith("/__init__.py") for name in sources),
            "tracked_script_paths": sum(name.startswith("scripts/") for name in tracked),
            "live_test_files": sum(name.startswith("tests/") and PurePosixPath(name).name.startswith("test_") for name in sources),
            "test_count_basis": "AST test functions/files only; no pytest collection or suite execution is performed by inventory",
            "open_decisions": [decision["id"] for decision in seed.get("decisions", []) if decision["status"] == "open"],
            "unresolved_role_paths": [entry["path"] for entry in entries if entry["role"] == "unowned_requires_decision"],
            "import_edge_count": len(graph["imports"]),
            "source_cycle_components": len(graph["source_strongly_connected_components"]),
            "dynamic_resolution_gap_count": len(graph["dynamic_resolution_gaps"]),
            "dynamic_gaps_by_context": {context: sum(gap["context"] == context for gap in graph["dynamic_resolution_gaps"]) for context in ("library", "script", "test", "historical_evidence")},
            "static_writer_gap_scope": "Conservative AST destinations only. These are not declarations that an artifact has no producer; exact producer/provenance review remains data/artifact_registry.json and W13.",
            "results_without_static_writer": [entry["path"] for entry in entries if entry.get("reproduction_status") == "producer_not_established_by_this_scan_consult_artifact_registry"],
            "retained_negative_claims": len(registry),
            "frozen_public_artifacts": len(frozen["records"]),
            "frozen_or_historical_hash_checks": len(frozen["integrity_checks"]),
            "verified_frozen_claim_assertions": sum(len(record["frozen_claim_checks"]) for record in registry),
            "protected_paths": len(protected),
            "completed_approved_cleanup_actions": len(completed),
            "pending_approved_cleanup_paths": sorted(set(cleanup["approved"]) - completed),
            "immutable_source_relocations": len(relocations),
            "architecture_complete": False,
            "inventory_consistency_is_not_scientific_or_release_acceptance": True,
        },
    }
    refreshed_workflow = deepcopy(workflow)
    for node in refreshed_workflow["nodes"]:
        added_tests = [item["path"] for item in workflow.get("reviewed_session_additions", []) if item["work_item"] == node["id"] and item["path"].startswith("tests/test_")]
        node["tests"] = list(dict.fromkeys(node["tests"] + added_tests))
    refreshed_workflow["topological_order"] = order
    refreshed_workflow["inventory_coverage"] = {node["id"]: sorted(entry["path"] for entry in entries if node["id"] in entry["work_items"]) for node in workflow["nodes"]}
    refreshed_workflow["dependency_edges"] = [{"prerequisite": predecessor, "deliverable": node["id"]} for node in workflow["nodes"] for predecessor in node.get("depends_on", [])]
    refreshed_workflow["open_decisions_by_work_item"] = {identity: [decision["id"] for decision in seed.get("decisions", []) if decision["status"] == "open" and identity in decision.get("blocks", [])] for identity in order}
    refreshed_deletions = validate_deletions(deletions, paths, protected, graph, frozen, negative, cleanup)
    return result, refreshed_workflow, refreshed_deletions


def check_snapshot(stored: dict, current: dict, stored_workflow: dict, current_workflow: dict, stored_deletions: dict, current_deletions: dict, verify_source_hashes=False) -> list[str]:
    errors = []
    stored_entries = stored.get("entries", [])
    names = [entry.get("path") for entry in stored_entries]
    if len(names) != len(set(names)):
        errors.append("Duplicate inventory paths")
    structural_fields = ("path", "tracked", "kind", "role", "purpose", "work_items", "accountable_owners", "baseline_git_object_id", "retained_negative_ids", "retirement_protected", "decision_ids", "source_identity")

    def project(entries):
        rows = []
        for entry in entries:
            row = {key: entry.get(key) for key in structural_fields}
            if "python" in entry:
                row["public_symbols"] = entry["python"]["declared_exports"]
                row["definition_names"] = [{key: definition[key] for key in ("name", "kind")} for definition in entry["python"]["definitions"]]
            rows.append(row)
        return rows

    def caller_signature(callers):
        return {
            "imports": semantic_graph(callers.get("static_imports", [])),
            "references": sorted({(item["source"], item["target"], item["kind"]) for item in callers.get("literal_or_document_references", [])}),
        }

    if stored.get("schema_version") != SCHEMA_VERSION or stored.get("kind") != "repository_inventory":
        errors.append("Invalid stored inventory schema/kind")
    if project(stored_entries) != project(current["entries"]):
        errors.append("Inventory coverage, roles, work ownership or retained-evidence membership changed; regenerate")
    if stored.get("scan", {}).get("semantic_import_graph_sha256") != current["scan"]["semantic_import_graph_sha256"]:
        errors.append("Static import topology changed; regenerate the use graph (ordinary body/hash/line edits do not trigger this check)")
    if stored.get("use_graph"):
        if digest(semantic_graph(stored["use_graph"]["imports"])) != stored["scan"].get("semantic_import_graph_sha256"):
            errors.append("Stored use graph does not match its semantic digest")
    if stored.get("summary") != current["summary"]:
        errors.append("Inventory counts, unresolved-role/reproduction lists or decision summary changed; regenerate")
    def portable_status(items):
        return [{key: item[key] for key in ("path", "tracked", "ignored", "contents_read")} for item in items]

    if portable_status(stored.get("special_root_artifacts", [])) != portable_status(current["special_root_artifacts"]):
        errors.append("Status-only root artifact Git metadata changed; regenerate")
    for field in ("baseline_commit", "baseline_tree_object_id", "tracked_pathset_sha256", "baseline_blob_catalog_sha256", "frozen_integrity_checks"):
        if stored.get("scan", {}).get(field) != current["scan"].get(field):
            errors.append(f"Scan authority or frozen integrity changed: {field}")
    if stored.get("retained_negative_registry") != current["retained_negative_registry"]:
        errors.append("Retained-negative registry expansion is stale")
    for field in ("source_identity_registry", "completed_cleanup_actions"):
        if stored.get(field) != current.get(field):
            errors.append(f"Verified historical identity/cleanup record changed: {field}")
    for field in ("nodes", "topological_order", "inventory_coverage", "dependency_edges", "open_decisions_by_work_item"):
        if stored_workflow.get(field) != current_workflow.get(field):
            errors.append(f"Workflow {field} is stale")
    for field in ("protected_exact_paths", "approval_request", "execution_verification"):
        if stored_deletions.get(field) != current_deletions.get(field):
            errors.append(f"Deletion approval/protection verification is stale: {field}")
    for old, new in zip(stored_deletions.get("proposals", []), current_deletions["proposals"]):
        if old.get("frozen_impact") != new.get("frozen_impact"):
            errors.append(f"Deletion frozen-impact review changed: {new['path']}")
        if caller_signature(old.get("callers", {})) != caller_signature(new["callers"]):
            errors.append(f"Deletion caller review changed: {new['path']}")
        if old.get("action_executed") is not new["action_executed"]:
            errors.append(f"Cleanup execution does not match its verified receipt: {new['path']}")
    if verify_source_hashes:
        expected = stored.get("scan", {}).get("source_hashes_sha256")
        if expected is None:
            errors.append("Exact source verification requires a generated --snapshot-hashes snapshot")
        elif expected != current["scan"].get("source_hashes_sha256"):
            actual = current["scan"].get("source_hashes_sha256", {})
            changed = sorted(name for name in set(expected) | set(actual) if expected.get(name) != actual.get(name))
            errors.append("Source bytes differ from the exact final snapshot: " + ", ".join(changed))
    return errors


def run(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Inventory tracked paths, actual static uses and protected evidence without importing ystwin or loading local secrets.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, help="Write/read the three inventory JSON documents here; default ROOT/data")
    parser.add_argument("--check", action="store_true", help="Read-only coverage, role, import-topology, DAG and frozen-integrity check; not scientific acceptance")
    parser.add_argument("--snapshot-hashes", action="store_true", help="Include exact public source hashes for a final external snapshot")
    parser.add_argument("--verify-source-hashes", action="store_true", help="With --check, also verify exact source bytes from --snapshot-hashes")
    parser.add_argument("--include-untracked", action="append", default=[], metavar="PATH", help="Explicitly inventory additional public source/CI paths without recursively reading untracked files")
    parser.add_argument("--require-resolved", action="store_true", help="Exit 2 while roles/decisions/dynamic local imports remain unresolved")
    arguments = parser.parse_args(argv)
    if arguments.verify_source_hashes and not arguments.check:
        parser.error("--verify-source-hashes requires --check")
    root = arguments.root.resolve()
    output = arguments.output_dir.resolve() if arguments.output_dir else root / "data"
    try:
        inputs = output if arguments.check else root / "data"
        documents = [parse_json(safe_file(inputs, name).read_bytes()) for name in DOCUMENTS]
        seed, workflow, deletions = documents
        include = arguments.include_untracked or seed.get("scan", {}).get("explicit_untracked_inclusions", [])
        current = build_inventory(root, seed, workflow, deletions, include, arguments.snapshot_hashes or arguments.verify_source_hashes)
        inventory, new_workflow, new_deletions = current
        if arguments.check:
            errors = check_snapshot(seed, inventory, workflow, new_workflow, deletions, new_deletions, arguments.verify_source_hashes)
            if errors:
                print(json.dumps({"check_passed": False, "errors": errors, "architecture_complete": False}, indent=2))
                return 1
        else:
            if not output.parent.is_dir():
                raise InventoryError("Output parent must already exist")
            output.mkdir(exist_ok=True)
            for name, document in zip(DOCUMENTS, current):
                target = output / name
                if target.is_symlink():
                    raise InventoryError(f"Refusing symlinked output: {name}")
                target.write_text(serialized(document) + "\n", encoding="utf-8")
        summary = deepcopy(inventory["summary"])
        summary["check_passed"] = True if arguments.check else None
        summary["mode"] = "read_only_check" if arguments.check else "generated"
        summary["exact_source_hashes_verified"] = bool(arguments.check and arguments.verify_source_hashes)
        print(json.dumps(summary, indent=2))
        if arguments.require_resolved and (summary["open_decisions"] or summary["unresolved_role_paths"] or summary["dynamic_resolution_gap_count"]):
            return 2
        return 0
    except (InventoryError, OSError, KeyError, TypeError) as exc:
        print(f"repository inventory refused: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(run())
