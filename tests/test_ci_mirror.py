"""The reviewer-facing CI mirror must match the live workflow definitions.

`.github/` cannot ship in the review archive (CODEOWNERS names the owner), so
`docs/CI_WORKFLOWS.md` mirrors the workflow YAML for reviewers. This test keeps
the mirror from silently drifting.
"""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MIRROR = ROOT / "docs" / "CI_WORKFLOWS.md"
WORKFLOWS = sorted((ROOT / ".github" / "workflows").glob("*.yml"))


def test_mirror_exists():
    assert MIRROR.is_file(), "docs/CI_WORKFLOWS.md is the only CI evidence reviewers get"


@pytest.mark.parametrize("workflow", WORKFLOWS, ids=lambda p: p.name)
def test_workflow_is_mirrored_verbatim(workflow):
    text = MIRROR.read_text(encoding="utf8")
    body = workflow.read_text(encoding="utf8").rstrip()
    assert f"`.github/workflows/{workflow.name}`" in text, "workflow not listed in the mirror"
    assert body in text, f"{workflow.name} differs from the mirrored copy; refresh the mirror"
