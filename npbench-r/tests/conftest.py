import textwrap
from pathlib import Path

import pytest
import yaml

from npbench.grading.config import load_grading_config


@pytest.fixture
def write(tmp_path):
    """Write a file under tmp_path, dedenting and stripping the leading newline."""

    def _write(relative: str, content: str) -> Path:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(content).lstrip("\n"))
        return path

    return _write


@pytest.fixture
def spec_from(tmp_path):
    """Build one ArtifactSpec by round-tripping a block through grading.yaml."""

    def _spec(block: dict):
        path = tmp_path / "grading.yaml"
        path.write_text(yaml.safe_dump({"version": 1, "artifacts": [block]}))
        return load_grading_config(path).artifacts[0]

    return _spec
