"""No script may pass a keyword argument the library no longer accepts.

`scripts/run_scenarios.py` passed ``maintenance_per_activity=`` to ``latent_constraints``.
That parameter was removed when the constant behind it was found to be dimensionally
incoherent -- it multiplied a promoter activity in RFU/OD/h and called the product
mmol ATP/gDW/h -- and the caller was never moved with it. **The script had been unable to
start since that refactor**, and its twelve tests passed the whole time, because both of
the ones that reach ``main()`` do so on the absent-model path and exit before the call.

A refactor that removes a parameter is a mechanical change with a mechanical check, and
this is it. It cannot catch a caller that passes the right names with the wrong meanings,
which is a separate problem for separate tests -- but it catches callers with resolvable
library identities that were simply left behind, which is what happened here.
"""

from __future__ import annotations

import ast
import builtins
from copy import copy
from dataclasses import dataclass
import importlib
import inspect
import pathlib
import symtable

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_UNBOUND = object()


def _accepted_keywords(obj) -> set[str] | None:
    """One resolved ystwin callable's keywords, or ``None`` when it cannot be checked.

    Imports and concrete receiver types establish identity; bare names do not. Inspect
    that callable alone, never merge namesakes. Only direct concrete annotations and
    return types are followed. Conflicting branches and loop/exception rebindings are
    unknown, as are container elements, union types and unannotated factories. Dynamic
    attributes, unavailable optional modules, expanded ``**kwargs`` and runtime mutation
    across call boundaries remain outside this gate's scope.
    A false negative here is a check that did not fire; a false positive would be a gate
    nobody could satisfy.
    """
    if not callable(obj) or (getattr(obj, "__module__", None) or "").partition(".")[0] != "ystwin":
        return None
    try:
        signature = inspect.signature(obj)
    except (ValueError, TypeError):
        return None
    if any(p.kind is p.VAR_KEYWORD for p in signature.parameters.values()):
        return None
    return {name for name, p in signature.parameters.items()
            if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)}


def _module(name):
    if not name or name.partition(".")[0] != "ystwin":
        return None
    try:
        return importlib.import_module(name)
    except Exception:  # an optional dependency this machine lacks is not this test's business
        return None


@dataclass(frozen=True)
class _Instance:
    cls: type


@dataclass(frozen=True)
class _Factory:
    returns: _Instance | None


def _annotation(annotation, env):
    if annotation is inspect.Signature.empty:
        return None
    if isinstance(annotation, str):
        try:
            annotation = ast.parse(annotation, mode="eval").body
        except SyntaxError:
            return None
    if isinstance(annotation, ast.Constant) and isinstance(annotation.value, str):
        return _annotation(annotation.value, env)
    obj = _resolve(annotation, env) if isinstance(annotation, ast.AST) else annotation
    return _Instance(obj) if isinstance(obj, type) else None


def _returns(obj):
    if isinstance(obj, _Factory):
        return copy(obj.returns)
    if isinstance(obj, type):
        return _Instance(obj)
    obj = inspect.unwrap(obj)
    if (inspect.isfunction(obj) or inspect.ismethod(obj)) and not (
            inspect.iscoroutinefunction(obj) or inspect.isgeneratorfunction(obj)):
        return _annotation(inspect.signature(obj).return_annotation, obj.__globals__)
    return None


def _attribute(obj, name):
    # METHODS TOO, and the omission was a false positive waiting to happen. A call
    # site reads `something.reaction_dg0(..., with_error=True)`; the old walk took the
    # bare attribute name, so a method sharing a name with a module-level function was
    # checked against the FUNCTION's signature and nothing else. On 2026-08-30
    # `bridge/equilibrator.py` hit exactly that -- a correct call to a correct method,
    # reported stale against a namesake it never touches. Merging signatures then hid
    # real errors and still mistook dict.update for ParticleFilter.update. Resolve the
    # receiver instead: each method belongs to its own class, not to a namesake's merge.
    instance = isinstance(obj, _Instance)
    owner = obj.cls if instance else obj
    if not (inspect.ismodule(owner) or isinstance(owner, type)):
        return None
    member = inspect.getattr_static(owner, name, None)
    if isinstance(member, (classmethod, staticmethod)):
        return member.__get__(None, owner)
    if instance and inspect.isfunction(member):
        return member.__get__(obj, owner)
    if instance and isinstance(member, property):
        return _returns(member.fget)
    if member is None and inspect.ismodule(owner):
        return _module(f"{owner.__name__}.{name}")
    return member


def _resolve(node, env):
    if isinstance(node, ast.Name):
        return env.get(node.id, vars(builtins).get(node.id))
    if isinstance(node, ast.Attribute):
        owner = _resolve(node.value, env)
        return env.get((id(owner), node.attr), _attribute(owner, node.attr))
    if isinstance(node, ast.Call):
        return _returns(_resolve(node.func, env))
    if isinstance(node, ast.Constant):
        return _Instance(type(node.value))
    cls = {ast.Dict: dict, ast.DictComp: dict, ast.List: list, ast.ListComp: list,
           ast.Tuple: tuple, ast.Set: set, ast.SetComp: set}.get(type(node))
    return _Instance(cls) if cls else None


def _same(left, right):
    if isinstance(left, _Instance) and isinstance(right, _Instance):
        return left.cls is right.cls
    if isinstance(left, _Factory) and isinstance(right, _Factory):
        return _same(left.returns, right.returns)
    if inspect.ismethod(left) and inspect.ismethod(right):
        return left.__func__ is right.__func__ and _same(left.__self__, right.__self__)
    return left is right


def _merge(*envs):
    return {key: envs[0].get(key) if all(_same(env.get(key), envs[0].get(key)) for env in envs)
            else None for key in set().union(*envs)}


class _Calls(ast.NodeVisitor):
    def __init__(self, path, text, stale):
        self.path, self.stale = path, stale
        try:
            parts = path.relative_to(_REPO / "src").with_suffix("").parts
        except ValueError:
            parts = ()
        dotted = ".".join(parts[:-1] if parts and parts[-1] == "__init__" else parts)
        self.package = dotted if path.name == "__init__.py" else dotted.rpartition(".")[0]
        self.owner = _module(dotted)
        self.env = dict(vars(self.owner)) if self.owner else {}
        self.outer = None
        self.tables = {}
        tables = [symtable.symtable(text, str(path), "exec")]
        while tables:
            table = tables.pop()
            self.tables[table.get_lineno(), table.get_name()] = table
            tables.extend(table.get_children())

    def scan(self, nodes):
        self.pending, self.changed = [], set()
        for node in nodes:
            self.visit(node)
        for node, owner, parameter_types in self.pending:
            child = copy(self)
            child.env = (self.outer if self.outer is not None else self.env).copy()
            for key in self.changed:
                child.env[key] = None
            table = self.tables[node.lineno, node.name]
            for symbol in table.get_symbols():
                if symbol.is_local() or symbol.is_assigned():
                    child.env[symbol.get_name()] = _UNBOUND
            child.env.update(parameter_types)
            positional = node.args.posonlyargs + node.args.args
            if isinstance(owner, type) and positional:
                member = inspect.getattr_static(owner, node.name, None)
                if isinstance(member, classmethod):
                    child.env[positional[0].arg] = owner
                elif inspect.isfunction(member):
                    child.env[positional[0].arg] = _Instance(owner)
            child.owner, child.outer = None, None
            child.scan(node.body)
        return self.env

    def branch(self, nodes, env=None):
        child = copy(self)
        child.env = (self.env if env is None else env).copy()
        result = child.scan(nodes)
        self.changed.update(child.changed)
        return result

    def bind(self, target, value):
        key = target.id if isinstance(target, ast.Name) else None
        if isinstance(target, ast.Attribute):
            owner = _resolve(target.value, self.env)
            if owner is None or owner is _UNBOUND:
                return
            key = (id(owner), target.attr)
        if key is not None:
            if self.env.get(key, _UNBOUND) is not _UNBOUND and not _same(self.env[key], value):
                self.changed.add(key)
            self.env[key] = value
        elif isinstance(target, (ast.Tuple, ast.List)):
            for item in target.elts:
                self.bind(item, None)

    def visit_Import(self, node):
        for alias in node.names:
            imported = _module(alias.name)
            name = alias.asname or alias.name.partition(".")[0]
            self.bind(ast.Name(id=name), imported if alias.asname else _module(name))

    def visit_ImportFrom(self, node):
        try:
            name = importlib.util.resolve_name("." * node.level + (node.module or ""), self.package)
        except (ImportError, ValueError):
            name = ""
        module = _module(name)
        for alias in node.names:
            if alias.name == "*":
                names = (getattr(module, "__all__", [n for n in vars(module) if not n.startswith("_")])
                         if module else list(self.env))
                for name in names:
                    self.bind(ast.Name(id=name), _attribute(module, name))
            else:
                self.bind(ast.Name(id=alias.asname or alias.name), _attribute(module, alias.name))

    def visit_Assign(self, node):
        self.visit(node.value)
        value = _resolve(node.value, self.env)
        for target in node.targets:
            self.bind(target, value)

    def visit_AnnAssign(self, node):
        if node.value:
            self.visit(node.value)
        value = _resolve(node.value, self.env) or _annotation(node.annotation, self.env)
        self.bind(node.target, value)

    def visit_NamedExpr(self, node):
        self.visit(node.value)
        self.bind(node.target, _resolve(node.value, self.env))

    def visit_AugAssign(self, node):
        self.generic_visit(node)
        self.bind(node.target, None)

    def visit_Delete(self, node):
        for target in node.targets:
            self.bind(target, None)

    def definition(self, node):
        obj = _attribute(self.owner, node.name)
        try:
            original = inspect.unwrap(obj)
            first_line = min([node.lineno] + [d.lineno for d in node.decorator_list])
            if (pathlib.Path(inspect.getsourcefile(original)).resolve() == self.path.resolve()
                    and inspect.getsourcelines(original)[1] == first_line):
                return obj
        except (TypeError, OSError):
            pass
        return None

    def visit_FunctionDef(self, node):
        for value in node.decorator_list + node.args.defaults + [v for v in node.args.kw_defaults if v]:
            self.visit(value)
        obj = self.definition(node)
        value = obj or (_Factory(_annotation(node.returns, self.env))
                        if not node.decorator_list and isinstance(node, ast.FunctionDef) else None)
        self.bind(ast.Name(id=node.name), value)
        parameter_types = {arg.arg: _annotation(arg.annotation, self.env)
                           for arg in node.args.posonlyargs + node.args.args + node.args.kwonlyargs}
        self.pending.append((node, self.owner, parameter_types))

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node):
        for value in node.decorator_list + node.bases + [kw.value for kw in node.keywords]:
            self.visit(value)
        obj = self.definition(node)
        self.bind(ast.Name(id=node.name), obj)
        child = copy(self)
        child.env, child.outer, child.owner = self.env.copy(), self.env, obj
        child.scan(node.body)

    def visit_If(self, node):
        self.visit(node.test)
        self.env.update(_merge(self.branch(node.body), self.branch(node.orelse)))

    def visit_IfExp(self, node):
        self.visit(node.test)
        self.env.update(_merge(self.branch([node.body]), self.branch([node.orelse])))

    def visit_BoolOp(self, node):
        self.visit(node.values[0])
        for value in node.values[1:]:
            self.env.update(_merge(self.env, self.branch([value])))

    def forget_writes(self, node):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            self.bind(ast.Name(id=node.name), None)
            return
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                self.bind(ast.Name(id=alias.asname or alias.name.partition(".")[0]), None)
        if isinstance(node, (ast.Name, ast.Attribute)) and isinstance(node.ctx, (ast.Store, ast.Del)):
            self.bind(node, None)
        if isinstance(node, (ast.ExceptHandler, ast.MatchAs, ast.MatchStar)) and node.name:
            self.bind(ast.Name(id=node.name), None)
        if isinstance(node, ast.MatchMapping) and node.rest:
            self.bind(ast.Name(id=node.rest), None)
        for child in ast.iter_child_nodes(node):
            self.forget_writes(child)

    def visit_For(self, node):
        self.visit(node.iter)
        self.forget_writes(node)
        self.branch(node.body)
        self.branch(node.orelse)

    visit_AsyncFor = visit_For

    def visit_While(self, node):
        self.forget_writes(node)
        self.visit(node.test)
        self.branch(node.body)
        self.branch(node.orelse)

    def visit_Try(self, node):
        self.forget_writes(node)
        body = self.branch(node.body)
        self.branch(node.orelse, body)
        for handler in node.handlers:
            if handler.type:
                self.visit(handler.type)
            self.branch(handler.body)
        for item in node.finalbody:
            self.visit(item)

    visit_TryStar = visit_Try

    def visit_Match(self, node):
        self.visit(node.subject)
        self.forget_writes(node)
        for case in node.cases:
            self.branch(([case.guard] if case.guard else []) + case.body)

    def visit_Lambda(self, node):
        for value in node.args.defaults + [v for v in node.args.kw_defaults if v]:
            self.visit(value)
        child = copy(self)
        child.env, child.changed = self.env.copy(), set()
        for arg in ast.walk(node.args):
            if isinstance(arg, ast.arg):
                child.bind(ast.Name(id=arg.arg), None)
        child.forget_writes(node.body)
        child.scan([node.body])

    def visit_ListComp(self, node):
        child = copy(self)
        child.env, child.changed = self.env.copy(), set()
        for generator in node.generators:
            child.visit(generator.iter)
            child.bind(generator.target, None)
            for condition in generator.ifs:
                child.visit(condition)
        for value in ([node.key, node.value] if isinstance(node, ast.DictComp) else [node.elt]):
            child.visit(value)

    visit_SetComp = visit_DictComp = visit_GeneratorExp = visit_ListComp

    def visit_With(self, node):
        for item in node.items:
            self.visit(item.context_expr)
            if item.optional_vars:
                self.bind(item.optional_vars, None)
        for item in node.body:
            self.visit(item)

    visit_AsyncWith = visit_With

    def visit_Call(self, node):
        obj = _resolve(node.func, self.env)
        accepted = _accepted_keywords(obj)
        if accepted is not None:
            for keyword in node.keywords:
                if keyword.arg and keyword.arg not in accepted:
                    try:
                        where = self.path.relative_to(_REPO)
                    except ValueError:      # a tmp_path fixture in the tests below
                        where = self.path
                    name = f"{obj.__module__}.{obj.__qualname__}"
                    self.stale.append(f"{where}:{node.lineno} {name}(..., {keyword.arg}=...) -- accepted: {sorted(accepted)}")
        self.generic_visit(node)


def _stale_keywords(files) -> list[str]:
    stale = []
    for path in files:
        try:
            text = path.read_text(encoding="utf-8")
            tree = ast.parse(text)
            _Calls(path, text, stale).scan(tree.body)
        except SyntaxError:
            continue
    return stale


class TestNoCallerWasLeftBehindByARefactor:
    def test_no_script_passes_a_keyword_the_library_does_not_accept(self):
        stale = _stale_keywords(sorted((_REPO / "scripts").rglob("*.py")))

        assert stale == []

    def test_no_library_module_passes_one_either(self):
        """The same check inward. A library that calls itself with a removed parameter is
        the same defect one layer down, and one layer harder to notice."""
        stale = _stale_keywords(sorted((_REPO / "src" / "ystwin").rglob("*.py")))

        assert stale == []


class TestTheCheckCanActuallyFail:
    """A gate that cannot fire is not a gate. These say the machinery notices."""

    def test_it_flags_a_keyword_that_is_not_in_the_signature(self, tmp_path):
        offender = tmp_path / "offender.py"
        offender.write_text("from ystwin.estimator import Observation\n"
                            "Observation(time_h=0, nonexistent_parameter_xyz=1)\n")

        stale = _stale_keywords([offender])
        assert len(stale) == 1
        assert "nonexistent_parameter_xyz=" in stale[0]

    def test_it_leaves_a_keyword_that_is_in_the_signature_alone(self, tmp_path):
        ok = tmp_path / "ok.py"
        ok.write_text("from ystwin.estimator import Observation\nObservation(time_h=0)\n")

        assert _stale_keywords([ok]) == []

    def test_a_name_the_library_does_not_export_is_not_policed(self, tmp_path):
        """The check only knows ystwin callables. A pandas or matplotlib call with a
        keyword it does not recognise must pass, or the gate becomes unsatisfiable."""
        other = tmp_path / "other.py"
        other.write_text("def write(frame, path):\n"
                         "    frame.to_csv(path, index=False, some_pandas_kwarg=1)\n")

        assert _stale_keywords([other]) == []


class TestCallableIdentity:
    def test_dictionary_update_is_not_particle_filter_update(self, tmp_path):
        source = tmp_path / "updates.py"
        source.write_text("from ystwin.estimator import ParticleFilter\n"
                          "def check(pf: ParticleFilter):\n"
                          "    protocol = {}\n"
                          "    protocol.update(run_complete=True)\n"
                          "    pf.update(run_complete=True)\n")

        stale = _stale_keywords([source])
        assert len(stale) == 1
        assert ":5 " in stale[0]
        assert "run_complete=" in stale[0]
        assert "observation" in stale[0]

    @pytest.mark.parametrize("source", [
        "from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "constrain(None, None, removed_parameter_xyz=1)\n",
        "import ystwin.bridge.latent_bridge\n"
        "ystwin.bridge.latent_bridge.latent_constraints(None, None, removed_parameter_xyz=1)\n",
        "import ystwin.bridge.latent_bridge as bridge\n"
        "bridge.latent_constraints(None, None, removed_parameter_xyz=1)\n",
        "from ystwin.bridge import latent_bridge as bridge\n"
        "bridge.latent_constraints(None, None, removed_parameter_xyz=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "constrain = latent_constraints\nconstrain(None, None, removed_parameter_xyz=1)\n",
        "from ystwin.estimator import ParticleFilter as PF\n"
        "def check(priors, optics):\n"
        "    pf = PF(priors, optics, 1.0, 0.0)\n"
        "    pf.update(removed_parameter_xyz=1)\n",
        "from ystwin.estimator import ParticleFilter as PF\n"
        "def check(pf: 'PF'):\n    pf.update(removed_parameter_xyz=1)\n",
        "from ystwin.estimator import ParticleFilter as PF\n"
        "def check(pf: PF):\n"
        "    update = pf.update\n    update(removed_parameter_xyz=1)\n",
        "from ystwin.estimator import ParticleFilter as PF\n"
        "def check(pf: PF):\n"
        "    pf.posterior().interval('biomass', removed_parameter_xyz=1)\n",
        "from ystwin.estimator import ParticleFilter as PF\n"
        "def make(priors, optics) -> PF:\n    return PF(priors, optics, 1.0, 0.0)\n"
        "def check(priors, optics):\n"
        "    make(priors, optics).update(removed_parameter_xyz=1)\n",
        "from ystwin.bridge.equilibrator import EquilibratorData\n"
        "EquilibratorData.load(removed_parameter_xyz=1)\n",
        "from ystwin.bridge.equilibrator import EquilibratorData\n"
        "thermo = EquilibratorData.load()\n"
        "thermo.reaction_dg0({}, removed_parameter_xyz=1)\n",
        "from ystwin.estimator import Observation as Record\nRecord(removed_parameter_xyz=1)\n",
        "import ystwin.estimator as estimator\n"
        "def check(pf: estimator.ParticleFilter):\n    pf.update(removed_parameter_xyz=1)\n",
        "from ystwin.mech.upr import UprChain\nUprChain.occupancy(None, removed_parameter_xyz=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "check = lambda: constrain(None, None, removed_parameter_xyz=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "[constrain(None, None, removed_parameter_xyz=1) for _ in range(1)]\n",
    ])
    def test_known_library_targets_are_still_checked(self, tmp_path, source):
        path = tmp_path / "library.py"
        path.write_text(source)

        stale = _stale_keywords([path])
        assert len(stale) == 1
        assert "removed_parameter_xyz=" in stale[0]

    @pytest.mark.parametrize("source", [
        "protocol = {}\nprotocol.update(run_complete=True)\n",
        "def check(receiver):\n    receiver.update(run_complete=True)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "def check(first, second):\n"
        "    first.update = constrain\n    second.update(run_complete=True)\n",
        "from collections import Counter\nCounter().update(run_complete=True)\n",
        "from types import SimpleNamespace as ParticleFilter\nParticleFilter(run_complete=True)\n",
        "def latent_constraints(**kwargs):\n    return kwargs\n"
        "latent_constraints(maintenance_per_activity=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "latent_constraints = dict\nlatent_constraints(maintenance_per_activity=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "def check(latent_constraints):\n    latent_constraints(maintenance_per_activity=1)\n",
        "from ystwin.estimator import ParticleFilter\n"
        "ParticleFilter = dict\nParticleFilter(run_complete=True)\n",
        "from types import SimpleNamespace\nimport ystwin.estimator as estimator\n"
        "estimator = SimpleNamespace(ParticleFilter=dict)\nestimator.ParticleFilter(run_complete=True)\n",
        "from ystwin.estimator import ParticleFilter\n"
        "def check(pf: ParticleFilter):\n    pf = {}\n    pf.update(run_complete=True)\n",
        "from ystwin.estimator import ParticleFilter\n"
        "def check(pf: ParticleFilter, use_dict):\n"
        "    if use_dict:\n        pf = {}\n    pf.update(run_complete=True)\n",
        "from ystwin.estimator import ParticleFilter\n"
        "def check(pf: ParticleFilter):\n"
        "    for pf in [{}]:\n        pf.update(run_complete=True)\n",
        "from ystwin.estimator import ParticleFilter\n"
        "def check(pf: ParticleFilter):\n"
        "    pf.update = dict\n    pf.update(run_complete=True)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "(lambda latent_constraints: latent_constraints(maintenance_per_activity=1))(dict)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "[latent_constraints(maintenance_per_activity=1) for latent_constraints in [dict]]\n",
        "from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "def check():\n    constrain(maintenance_per_activity=1)\nconstrain = dict\ncheck()\n",
        "constrain = dict\ndef check():\n    constrain(maintenance_per_activity=1)\ncheck()\n"
        "from ystwin.bridge.latent_bridge import latent_constraints as constrain\n",
        "def check(fail):\n    constrain = dict\n    try:\n        fail()\n"
        "    except Exception:\n        from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "    constrain(maintenance_per_activity=1)\n",
        "def check(again):\n    constrain = dict\n    while again():\n"
        "        from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
        "    constrain(maintenance_per_activity=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "def check(flag):\n    constrain = dict\n"
        "    flag and (constrain := latent_constraints)\n"
        "    constrain(maintenance_per_activity=1)\n",
        "from ystwin.bridge.latent_bridge import latent_constraints\n"
        "def check(flag):\n    constrain = dict\n"
        "    (constrain := latent_constraints) if flag else None\n"
        "    constrain(maintenance_per_activity=1)\n",
    ])
    def test_non_library_or_unresolved_targets_are_not_guessed(self, tmp_path, source):
        path = tmp_path / "other.py"
        path.write_text(source)

        assert _stale_keywords([path]) == []

    def test_rebinding_is_local_to_its_scope(self, tmp_path):
        source = tmp_path / "scopes.py"
        source.write_text("from ystwin.bridge.latent_bridge import latent_constraints as constrain\n"
                          "def local():\n"
                          "    constrain = dict\n    constrain(maintenance_per_activity=1)\n"
                          "constrain(None, None, maintenance_per_activity=1)\n")

        stale = _stale_keywords([source])
        assert len(stale) == 1
        assert ":5 " in stale[0]
        assert "maintenance_per_activity=" in stale[0]

    def test_namesake_functions_and_methods_keep_distinct_signatures(self, tmp_path):
        source = tmp_path / "namesakes.py"
        source.write_text("from ystwin.bridge.equilibrator import EquilibratorData\n"
                          "from ystwin.bridge.thermodynamic import reaction_dg0\n"
                          "def check(thermo: EquilibratorData):\n"
                          "    thermo.reaction_dg0({}, with_error=True)\n"
                          "    reaction_dg0({}, None, with_error=True)\n")

        stale = _stale_keywords([source])
        assert len(stale) == 1
        assert ":5 " in stale[0]
        assert "with_error=" in stale[0]

    def test_kwargs_on_one_classmethod_do_not_exempt_its_namesake(self, tmp_path):
        source = tmp_path / "classmethods.py"
        source.write_text("from ystwin.bridge.equilibrator import EquilibratorData, PreferredEnergies\n"
                          "from ystwin.bridge.thermodynamic import ThermodynamicData\n"
                          "PreferredEnergies.load(thermodb={})\n"
                          "ThermodynamicData.load(thermodb={})\n"
                          "EquilibratorData.load(thermodb={})\n")

        stale = _stale_keywords([source])
        assert len(stale) == 1
        assert ":5 " in stale[0]
        assert "thermodb=" in stale[0]

    def test_a_stale_internal_self_call_is_detected(self):
        path = _REPO / "src" / "ystwin" / "estimator.py"
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text)
        call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call)
                    and ast.unparse(node.func) == "self.update")
        call.keywords.append(ast.keyword(arg="removed_parameter_xyz", value=ast.Constant(1)))
        stale = []

        _Calls(path, text, stale).scan(tree.body)

        assert len(stale) == 1
        assert "ystwin.estimator.ParticleFilter.update" in stale[0]
        assert "removed_parameter_xyz=" in stale[0]

    @pytest.mark.parametrize("filename", ["equilibrator.py", "__init__.py"])
    def test_relative_imports_use_the_scanned_package(self, tmp_path, monkeypatch, filename):
        source = tmp_path / "src" / "ystwin" / "bridge" / filename
        source.parent.mkdir(parents=True)
        source.write_text("from .latent_bridge import latent_constraints as constrain\n"
                          "from ..estimator import ParticleFilter as PF\n"
                          "constrain(None, None, maintenance_per_activity=1)\n"
                          "def check(pf: PF):\n    pf.update(run_complete=True)\n")
        monkeypatch.setattr(__import__(__name__), "_REPO", tmp_path)

        stale = _stale_keywords([source])
        assert len(stale) == 2
        assert any(":3 " in item and "maintenance_per_activity=" in item for item in stale)
        assert any(":5 " in item and "run_complete=" in item for item in stale)
