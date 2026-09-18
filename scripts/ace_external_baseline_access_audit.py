#!/usr/bin/env python3
"""Record public-baseline access status for the funding-grade ACE upgrade."""

from __future__ import annotations

import json
import subprocess
import urllib.error
import urllib.request
import zipfile
from datetime import date
from pathlib import Path
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
P5_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "05-peptide-ai-benchmark"
AHTPIN_DIR = REPO_ROOT / "data" / "external-baselines" / "ahtpin"
AHTPIN_ZIP = AHTPIN_DIR / "AHTpin Dataset-20260508T063233Z-3-001.zip"
OUT_JSON = P5_DIR / "external-baseline-access-audit.json"
OUT_REPORT = P5_DIR / "external-baseline-access-audit.md"

PUBLIC_BASELINES = (
    {
        "baseline_id": "peptideranker_official_server",
        "role": "general bioactive-peptide predictor baseline",
        "primary_url": "https://bioware.ucd.ie/~compass/biowareweb/Server_pages/peptideranker.php",
        "documentation_url": "https://bioware.ucd.ie/~compass/biowareweb/Server_pages/help/peptideranker/help.php",
        "citation_url": "https://pubmed.ncbi.nlm.nih.gov/23056189/",
        "paper_or_tool_claim": "general bioactive-peptide ranking, not ACE-specific activity prediction",
    },
    {
        "baseline_id": "ahtpin_official_server_and_github_snapshot",
        "role": "ACE/antihypertensive peptide public predictor baseline candidate",
        "primary_url": "https://webs.iiitd.edu.in/raghava/ahtpin/",
        "documentation_url": "https://github.com/raghavagps/AHTpin",
        "citation_url": "https://doi.org/10.1038/srep12512",
        "paper_or_tool_claim": "antihypertensive peptide prediction/screening/design platform",
    },
    {
        "baseline_id": "uwm_bioactivity_predictor_index",
        "role": "baseline discovery index",
        "primary_url": "https://biochemia.uwm.edu.pl/en/bioactivity-prediction-2/",
        "documentation_url": "https://biochemia.uwm.edu.pl/en/bioactivity-prediction-2/",
        "citation_url": "",
        "paper_or_tool_claim": "curated links to peptide bioactivity predictors and docking tools",
    },
)


def request_status(url: str, timeout: int = 15) -> dict[str, Any]:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "peptide-lab-audit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return {
                "url": url,
                "status": "reachable",
                "http_status": response.status,
                "content_type": response.headers.get("content-type", ""),
            }
    except urllib.error.HTTPError as exc:
        return {
            "url": url,
            "status": "http_error",
            "http_status": exc.code,
            "error": str(exc),
        }
    except Exception as exc:  # noqa: BLE001 - audit should record access failures, not crash.
        return {
            "url": url,
            "status": "access_failed",
            "http_status": "",
            "error": f"{type(exc).__name__}: {exc}",
        }


def git_value(args: Sequence[str], cwd: Path) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=cwd, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return ""


def ahtpin_snapshot() -> dict[str, Any]:
    files = sorted(str(path.relative_to(AHTPIN_DIR)) for path in AHTPIN_DIR.rglob("*") if path.is_file()) if AHTPIN_DIR.exists() else []
    visible_files = [path for path in files if not path.startswith(".git/")]
    zip_members: list[str] = []
    if AHTPIN_ZIP.exists():
        with zipfile.ZipFile(AHTPIN_ZIP) as archive:
            zip_members = sorted(archive.namelist())
    executable_like = [
        path
        for path in visible_files
        if Path(path).suffix.lower() in {".py", ".r", ".pl", ".sh", ".jar", ".model", ".pkl", ".pickle"}
    ]
    return {
        "local_path": str(AHTPIN_DIR),
        "repo_exists": AHTPIN_DIR.exists(),
        "remote": git_value(["remote", "get-url", "origin"], AHTPIN_DIR) if AHTPIN_DIR.exists() else "",
        "commit": git_value(["rev-parse", "HEAD"], AHTPIN_DIR) if AHTPIN_DIR.exists() else "",
        "visible_files": visible_files,
        "zip_members": zip_members,
        "executable_or_model_files": executable_like,
        "local_runnable_predictor_found": bool(executable_like),
    }


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def status_text(entry: Mapping[str, Any]) -> str:
    access = entry.get("primary_access", {})
    status = access.get("status", "")
    http_status = access.get("http_status", "")
    return f"{status}" + (f" ({http_status})" if http_status else "")


def report_text(payload: Mapping[str, Any]) -> str:
    ahtpin = payload["ahtpin_local_snapshot"]
    rows = []
    for entry in payload["baselines"]:
        runnable = entry.get("runnable_status", "")
        rows.append(
            [
                entry["baseline_id"],
                entry["role"],
                status_text(entry),
                runnable,
            ]
        )
    table = markdown_table(["Baseline", "Role", "Access", "Runnable for v6 now"], rows)
    return f"""# External Baseline Access Audit

Date: {payload['date']}

## Purpose

Separate real public-baseline evidence from baseline names that were merely identified. This prevents the ACE v6 funding story from leaning on unavailable or unrun comparators.

## Access Results

{table}

## AHTpin Local Snapshot

- Remote: `{ahtpin['remote']}`
- Commit: `{ahtpin['commit']}`
- Local files excluding `.git`: `{len(ahtpin['visible_files'])}`
- Dataset zip members: `{len(ahtpin['zip_members'])}`
- Executable/model files found: `{len(ahtpin['executable_or_model_files'])}`

The cloned AHTpin snapshot contains documentation and a dataset archive, but no local trained predictor or executable scoring pipeline. It is therefore a baseline candidate, not a completed scored baseline.

## Claim Boundary

- Current manuscript-ready comparator remains the transparent local PeptideRanker-style composition/length baseline.
- Official PeptideRanker and AHTpin scores must not be claimed until those tools are actually run on the frozen v6 universe.
- AHTpin is the stronger target-relevant baseline to pursue first because it is antihypertensive/ACE-oriented; PeptideRanker is useful as a general bioactivity comparator.

## Next Baseline Work

1. Obtain a batch route or author-supported scoring route for AHTpin.
2. Run the unchanged v6 blinded candidate universe through AHTpin and compare top-1% recovery plus matched-decoy AUROC.
3. Retry official PeptideRanker only as a general bioactivity baseline, not as an ACE-specific comparator.
4. Keep the local PeptideRanker-style baseline in the paper as a reproducible fallback and clearly label it as local.
"""


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def main() -> int:
    baselines: list[dict[str, Any]] = []
    for baseline in PUBLIC_BASELINES:
        entry = dict(baseline)
        entry["primary_access"] = request_status(baseline["primary_url"])
        entry["documentation_access"] = request_status(baseline["documentation_url"])
        if baseline["baseline_id"] == "ahtpin_official_server_and_github_snapshot":
            snapshot = ahtpin_snapshot()
            entry["runnable_status"] = (
                "yes" if snapshot["local_runnable_predictor_found"] else "no_local_predictor_or_model_in_snapshot"
            )
        elif baseline["baseline_id"] == "peptideranker_official_server":
            entry["runnable_status"] = "not_run_from_this_host"
        else:
            entry["runnable_status"] = "index_only"
        baselines.append(entry)

    payload = {
        "date": date.today().isoformat(),
        "baselines": baselines,
        "ahtpin_local_snapshot": ahtpin_snapshot(),
    }
    write_json(OUT_JSON, payload)
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.write_text(report_text(payload), encoding="utf-8")
    print(json.dumps({"baselines": len(baselines), "report": str(OUT_REPORT)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
