#!/usr/bin/env python3
"""Write a reproducibility manifest for ACE submission artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import PAPER_DIR, save_json  # noqa: E402

DEFAULT_OUT_JSON = PAPER_DIR / "submission" / "reproducibility-manifest.json"
DEFAULT_OUT_MD = PAPER_DIR / "submission" / "reproducibility-manifest.md"

ARTIFACT_GLOBS = (
    "known-positive-*.csv",
    "frozen-recipes/**/*.json",
    "tables/*.csv",
    "tables/*.json",
    "figures/submission/*.svg",
    "submission/*.md",
    "experiments/2026-08-29-v4-independent-validation-001/frozen-rerank/*/frozen-rerank-summary.csv",
    "experiments/2026-08-29-v5-failure-mode-validation-001/frozen-rerank/*/frozen-rerank-summary.csv",
    "experiments/2026-08-29-v6-independent-validation-001/frozen-rerank/*/frozen-rerank-summary.csv",
    "experiments/2026-08-29-v6-independent-validation-001/v6-*.csv",
    "experiments/2026-08-29-v6-independent-validation-001/v6-*.json",
    "experiments/2026-08-29-v6-independent-validation-001/v6-*.md",
    "experiments/2026-08-29-v*-*/external-baselines/**/*.csv",
    "experiments/2026-08-29-v*-*/external-baselines/**/*.json",
    "experiments/2026-08-29-v*-*/external-baselines/**/*.md",
)
SCRIPT_PATHS = (
    REPO_ROOT / "scripts" / "ace_build_blinded_universe.py",
    REPO_ROOT / "scripts" / "ace_audit_panel_overlap.py",
    REPO_ROOT / "scripts" / "ace_apply_frozen_rerank_recipe.py",
    REPO_ROOT / "scripts" / "ace_run_scoring_mode_sweep.py",
    REPO_ROOT / "scripts" / "ace_train_rank_fusion_model.py",
    REPO_ROOT / "scripts" / "ace_build_structured_assay_metadata.py",
    REPO_ROOT / "scripts" / "ace_run_peptideranker_style_baseline.py",
    REPO_ROOT / "scripts" / "ace_motif_neighborhood_audit.py",
    REPO_ROOT / "scripts" / "ace_build_submission_artifacts.py",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", default=str(DEFAULT_OUT_JSON))
    parser.add_argument("--out-md", default=str(DEFAULT_OUT_MD))
    return parser.parse_args()


def command_output(args: list[str]) -> str:
    return subprocess.check_output(args, cwd=REPO_ROOT, text=True).strip()


def git_metadata() -> dict[str, object]:
    head = command_output(["git", "rev-parse", "HEAD"])
    status = command_output(["git", "status", "--short"])
    branch = command_output(["git", "rev-parse", "--abbrev-ref", "HEAD"])
    return {
        "head_commit": head,
        "branch": branch,
        "dirty": bool(status),
        "status_short": status.splitlines(),
        "archive_commit_status": "pending_clean_commit" if status else "fixed_clean_commit",
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_paths() -> list[Path]:
    paths: set[Path] = set()
    for pattern in ARTIFACT_GLOBS:
        paths.update(path for path in PAPER_DIR.glob(pattern) if path.is_file())
    paths.update(path for path in SCRIPT_PATHS if path.exists())
    return sorted(paths)


def relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def markdown(manifest: dict[str, object]) -> str:
    git = manifest["git"]
    artifacts = manifest["artifacts"]
    assert isinstance(git, dict)
    assert isinstance(artifacts, list)
    lines = [
        "# ACE Reproducibility Manifest",
        "",
        f"Date: {manifest['date']}",
        "",
        "## Git State",
        "",
        f"- Branch: `{git['branch']}`",
        f"- HEAD: `{git['head_commit']}`",
        f"- Dirty worktree: `{git['dirty']}`",
        f"- Archive status: `{git['archive_commit_status']}`",
        "",
        "A journal archive should be cut from a clean commit. This manifest records current artifact hashes and explicitly marks whether the worktree is clean.",
        "",
        "## Artifact Hashes",
        "",
        "| Path | SHA256 |",
        "| --- | --- |",
    ]
    for item in artifacts:
        if isinstance(item, dict):
            lines.append(f"| `{item['path']}` | `{item['sha256']}` |")
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    artifacts = [
        {
            "path": relative(path),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        }
        for path in artifact_paths()
    ]
    manifest = {
        "date": date.today().isoformat(),
        "git": git_metadata(),
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
        "claim_boundary": "This manifest is a pre-submission reproducibility record. Final archival release should use a clean commit and DOI.",
    }
    save_json(Path(args.out_json), manifest)
    Path(args.out_md).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_md).write_text(markdown(manifest), encoding="utf-8")
    print(args.out_json)
    print(args.out_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
