"""Every refusal a build can make must be demonstrated by a test that names it.

A template's `_assert_has_content` is the only thing standing between a campaign
whose rungs measure nothing and a campaign that ships. Each of its refusals is a
claim about what the data must look like, and a claim no test exercises is a
claim nobody has checked -- the same defect this project keeps finding in other
people's declarations, turned on its own.

Two failure modes matter and they are different:

- A refusal is ADDED without a test. Nothing notices: the build still passes on
  real gold, the oracle still scores 1.0, and the refusal might be unreachable or
  inverted. This module is the check for that.
- A refusal's condition becomes UNREACHABLE while its test still exists. The test
  notices, because `pytest.raises` fails with "DID NOT RAISE". That happened once
  in this project, to a strata check whose guard counted a census key that was
  always present, so the two mechanisms together cover both directions.

The check is static and order-independent on purpose. Tracing the suite and
asserting every `raise` line executed would be more precise, but it would depend
on the whole suite running in one process and on test order, and this project
runs its tests under `pytest-randomly`.

What it enforces: for every `raise` inside `_assert_has_content`, the template's
campaign test file carries a pattern that matches that refusal's message -- either
in a `pytest.raises(..., match="...")` directly, or as a literal in a
`@pytest.mark.parametrize` list feeding a `match=` variable, which is how the
parametrised refusal tests here are written. That is why every refusal message in this project
carries a distinctive phrase -- the phrase is the contract between the refusal
and its test.
"""

from __future__ import annotations

import ast
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "src" / "npbench_c" / "templates"
UNIT = ROOT / "tests" / "unit"

#: Which test file carries each template's refusal tests. Declared rather than
#: guessed, so a template whose tests move has to say so here and a template with
#: no entry fails the discovery test below rather than being skipped silently.
TEST_FILES = {
    "annotation_transfer": "test_campaign_transfer.py",
    "bgc_detection": "test_campaign_bgcdetect.py",
    "chembl_selectivity": "test_campaign_selectivity.py",
    "chemical_space": "test_campaign_chemspace.py",
    "domain_architecture": "test_campaign_arch.py",
    "ec_domain_audit": "test_campaign_ecaudit.py",
    "frameshift_response": "test_campaign_frameshift.py",
    "gcf_cutoff": "test_campaign_gcfcutoff.py",
    "gene_calling": "test_campaign_genecall.py",
    "kinetics_consistency": "test_campaign_kinetics.py",
    "ms2_dereplication": "test_campaign_derep.py",
    "remote_homology": "test_campaign_homology.py",
    "ripp_precursor": "test_campaign_ripp.py",
}


def _literal_parts(node: ast.AST) -> list[str]:
    """The literal text of a message expression, with holes for interpolations.

    An f-string's `{...}` parts carry no fixed text, so they become a space: a
    test's `match` pattern can only ever anchor on the literal spans around them.
    """
    parts: list[str] = []

    def walk(n: ast.AST) -> None:
        if isinstance(n, ast.Constant) and isinstance(n.value, str):
            parts.append(n.value)
        elif isinstance(n, ast.JoinedStr):
            for value in n.values:
                walk(value)
        elif isinstance(n, ast.BinOp):
            walk(n.left)
            walk(n.right)
        elif isinstance(n, ast.FormattedValue):
            parts.append(" ")

    walk(node)
    return parts


def _assert_function(path: pathlib.Path) -> ast.FunctionDef | None:
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.FunctionDef) and node.name == "_assert_has_content":
            return node
    return None


def _name_literals(function: ast.FunctionDef, name: str) -> list[str]:
    """Literal strings a local name can hold inside this function.

    Several templates loop over `(field, complaint)` pairs and raise the
    complaint, so the message is a Name rather than a literal at the raise site.
    The strings are still right there in the function, and a refusal written that
    way is no less a refusal, so they are collected rather than excused.
    """
    out: list[str] = []
    for node in ast.walk(function):
        if isinstance(node, ast.Assign):
            targets, source = list(node.targets), node.value
        elif isinstance(node, (ast.For, ast.comprehension)):
            targets, source = [node.target], node.iter
        else:
            continue
        for target in targets:
            if isinstance(target, ast.Name):
                if target.id != name:
                    continue
                for literal in ast.walk(source):
                    if (isinstance(literal, ast.Constant)
                            and isinstance(literal.value, str)):
                        out.append(literal.value)
            elif isinstance(target, (ast.Tuple, ast.List)):
                names = [e.id if isinstance(e, ast.Name) else None
                         for e in target.elts]
                if name not in names:
                    continue
                # Only the tuple POSITION this name unpacks from. Collecting
                # every literal in the iterable would demand a test for the
                # field names sitting alongside the messages, which are not
                # refusals and have nothing to name.
                index = names.index(name)
                for row in ast.walk(source):
                    if not isinstance(row, (ast.Tuple, ast.List)):
                        continue
                    if len(row.elts) != len(names):
                        continue
                    out.extend(_literal_parts(row.elts[index]))
    return out


def refusals(template: str) -> list[tuple[int, str]]:
    """(line, flattened message) for every raise inside `_assert_has_content`."""
    function = _assert_function(TEMPLATES / template / "build.py")
    if function is None:
        return []
    out: list[tuple[int, str]] = []
    for node in ast.walk(function):
        if not isinstance(node, ast.Raise) or node.exc is None:
            continue
        call = node.exc
        if not isinstance(call, ast.Call) or not call.args:
            continue
        argument = call.args[0]
        if isinstance(argument, ast.Name):
            candidates = _name_literals(function, argument.id)
        else:
            candidates = _literal_parts(argument)
        # One raise, one message: the parts are joined for a literal message and
        # kept apart for a name that can hold several, so a test naming any of
        # them covers the raise.
        if isinstance(argument, ast.Name):
            out.extend((node.lineno, " ".join(c.split())) for c in candidates)
        else:
            out.append((node.lineno, " ".join("".join(candidates).split())))
    return out


def _is_parametrize(node: ast.AST) -> bool:
    """`@pytest.mark.parametrize(...)`, however it is spelled."""
    if not isinstance(node, ast.Call):
        return False
    func = node.func
    while isinstance(func, ast.Attribute):
        if func.attr == "parametrize":
            return True
        func = func.value
    return isinstance(func, ast.Name) and func.id == "parametrize"


def _is_raises(node: ast.AST) -> bool:
    if not isinstance(node, ast.Call):
        return False
    func = node.func
    if isinstance(func, ast.Attribute):
        return func.attr == "raises"
    return isinstance(func, ast.Name) and func.id == "raises"


def _declared_patterns(template: str) -> list[str]:
    """Every string a refusal test could be matching on, from both idioms.

    `pytest.raises(Err, match="...")` is the direct one. The other is a
    parametrised refusal test -- `@pytest.mark.parametrize([(mutate, message)])`
    with `match=message` -- where the pattern is a literal in the decorator and a
    variable at the call site. Both are in use here, and a regex over `match=`
    sees only the first, which is why this reads the syntax instead.
    """
    tree = ast.parse((UNIT / TEST_FILES[template]).read_text())
    out: list[str] = []
    for node in ast.walk(tree):
        if _is_raises(node):
            for keyword in node.keywords:
                if keyword.arg == "match":
                    out.extend(_literal_parts(keyword.value))
        elif _is_parametrize(node):
            for argument in node.args:
                for literal in ast.walk(argument):
                    if (isinstance(literal, ast.Constant)
                            and isinstance(literal.value, str)):
                        out.append(literal.value)
    return [pattern for pattern in out if pattern.strip()]


ALL = sorted(TEST_FILES)


def test_every_template_with_refusals_declares_a_test_file():
    """A template that refuses and is not listed above would be skipped by every
    check here, so the listing has to be complete rather than best-effort."""
    missing = []
    for build in sorted(TEMPLATES.glob("*/build.py")):
        template = build.parent.name
        if _assert_function(build) is None:
            continue
        if template not in TEST_FILES:
            missing.append(template)
    assert not missing, (
        f"these templates refuse but name no test file in TEST_FILES: {missing}")


@pytest.mark.parametrize("template", ALL)
def test_the_declared_test_file_exists(template):
    assert (UNIT / TEST_FILES[template]).is_file()


@pytest.mark.parametrize("template", ALL)
def test_the_template_actually_refuses(template):
    """Guards the guard: if a template's `_assert_has_content` were emptied, the
    coverage check below would pass vacuously."""
    assert refusals(template), template


@pytest.mark.parametrize("template", ALL)
def test_every_refusal_is_named_by_a_test(template):
    patterns = _declared_patterns(template)
    assert patterns, f"{TEST_FILES[template]} has no pytest.raises(match=...) at all"
    unnamed = []
    for lineno, message in refusals(template):
        if not message:
            unnamed.append((lineno, "<no literal text in the message>"))
            continue
        if not any(re.search(pattern, message) for pattern in patterns):
            unnamed.append((lineno, message[:88]))
    assert not unnamed, (
        f"{template}: {len(unnamed)} refusal(s) in _assert_has_content that no "
        f"pytest.raises(match=...) in {TEST_FILES[template]} names:\n"
        + "\n".join(f"  build.py:{line}  {msg}" for line, msg in unnamed))
