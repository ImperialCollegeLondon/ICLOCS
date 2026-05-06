"""Update persisted GitHub usage metrics for README badges.

This script intentionally lives outside the MATLAB source tree. GitHub only
exposes clone/view traffic for the previous 14 days, so the workflow stores
the per-day values it sees and accumulates from the first successful run.
"""

from __future__ import annotations

import json
import os
import re
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


API_ROOT = "https://api.github.com"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def infer_repository() -> str | None:
    config_path = Path(".git/config")
    if not config_path.exists():
        return None

    config = config_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"url\s*=\s*(?:git@github\.com:|https://github\.com/)([^/\s]+/[^.\s]+)(?:\.git)?", config)
    return match.group(1) if match else None


def request_json(path: str, token: str | None) -> tuple[Any, dict[str, str]]:
    request = urllib.request.Request(f"{API_ROOT}{path}")
    request.add_header("Accept", "application/vnd.github+json")
    request.add_header("X-GitHub-Api-Version", "2022-11-28")
    if token:
        request.add_header("Authorization", f"Bearer {token}")

    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            headers = {key.lower(): value for key, value in response.headers.items()}
            return json.load(response), headers
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code}: {detail[:500]}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error: {exc.reason}") from exc


def next_page_path(headers: dict[str, str]) -> str | None:
    link_header = headers.get("link", "")
    for part in link_header.split(","):
        if 'rel="next"' not in part:
            continue
        match = re.search(r"<https://api\.github\.com([^>]+)>", part)
        if match:
            return match.group(1)
    return None


def request_all_pages(path: str, token: str | None) -> list[Any]:
    items: list[Any] = []
    next_path: str | None = path
    while next_path:
        page, headers = request_json(next_path, token)
        if not isinstance(page, list):
            raise RuntimeError(f"Expected a list response for {next_path}")
        items.extend(page)
        next_path = next_page_path(headers)
    return items


def default_metrics(repository: str) -> dict[str, Any]:
    return {
        "repository": repository,
        "last_updated": None,
        "source": "GitHub REST API",
        "limitations": [
            "Clone and view traffic are only available from GitHub for the previous 14 days.",
            "Recorded clone totals start from the first successful workflow run.",
            "Release download totals count GitHub release assets, not source-code ZIP downloads.",
        ],
        "clones": {
            "total": 0,
            "unique_daily_total": 0,
            "days": {},
            "latest_14_days": {"count": 0, "uniques": 0},
        },
        "views": {
            "total": 0,
            "unique_daily_total": 0,
            "days": {},
            "latest_14_days": {"count": 0, "uniques": 0},
        },
        "release_asset_downloads": {
            "total": 0,
            "assets": [],
        },
        "warnings": [],
    }


def load_metrics(path: Path, repository: str) -> dict[str, Any]:
    if not path.exists():
        return default_metrics(repository)

    with path.open(encoding="utf-8") as handle:
        metrics = json.load(handle)

    baseline = default_metrics(repository)
    for key, value in baseline.items():
        metrics.setdefault(key, value)
    metrics["repository"] = repository
    return metrics


def update_daily_metric(metrics: dict[str, Any], metric_name: str, payload: dict[str, Any], item_key: str) -> None:
    metric = metrics.setdefault(metric_name, {})
    days = metric.setdefault("days", {})
    for item in payload.get(item_key, []):
        date = str(item.get("timestamp", ""))[:10]
        if not date:
            continue
        days[date] = {
            "count": int(item.get("count", 0)),
            "uniques": int(item.get("uniques", 0)),
        }

    ordered_days = dict(sorted(days.items()))
    metric["days"] = ordered_days
    metric["total"] = sum(day["count"] for day in ordered_days.values())
    metric["unique_daily_total"] = sum(day["uniques"] for day in ordered_days.values())
    metric["latest_14_days"] = {
        "count": int(payload.get("count", 0)),
        "uniques": int(payload.get("uniques", 0)),
    }


def update_release_downloads(metrics: dict[str, Any], repository: str, token: str | None) -> None:
    releases = request_all_pages(f"/repos/{repository}/releases?per_page=100", token)
    assets: list[dict[str, Any]] = []

    for release in releases:
        tag_name = release.get("tag_name") or release.get("name") or "untagged"
        for asset in release.get("assets", []):
            assets.append(
                {
                    "release": tag_name,
                    "name": asset.get("name"),
                    "download_count": int(asset.get("download_count", 0)),
                    "browser_download_url": asset.get("browser_download_url"),
                    "updated_at": asset.get("updated_at"),
                }
            )

    assets.sort(key=lambda asset: (asset["release"], asset["name"] or ""))
    metrics["release_asset_downloads"] = {
        "total": sum(asset["download_count"] for asset in assets),
        "assets": assets,
    }


def compact_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M".replace(".0M", "M")
    if value >= 1_000:
        return f"{value / 1_000:.1f}k".replace(".0k", "k")
    return str(value)


def badge(label: str, value: int, color: str) -> dict[str, Any]:
    return {
        "schemaVersion": 1,
        "label": label,
        "message": compact_count(value),
        "color": color,
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    repository = os.environ.get("GITHUB_REPOSITORY") or os.environ.get("REPOSITORY") or infer_repository()
    if not repository or "/" not in repository:
        print("Unable to determine GitHub repository. Set GITHUB_REPOSITORY=owner/name.", file=sys.stderr)
        return 2

    metrics_path = Path(os.environ.get("METRICS_PATH", "docs/metrics/github-usage.json"))
    clones_badge_path = Path(os.environ.get("CLONES_BADGE_PATH", "docs/metrics/clones-badge.json"))
    releases_badge_path = Path(os.environ.get("RELEASES_BADGE_PATH", "docs/metrics/release-downloads-badge.json"))
    token = os.environ.get("TRAFFIC_TOKEN") or os.environ.get("GITHUB_TOKEN")

    metrics = load_metrics(metrics_path, repository)
    warnings: list[str] = []

    if token:
        for metric_name, item_key in (("clones", "clones"), ("views", "views")):
            try:
                payload, _ = request_json(f"/repos/{repository}/traffic/{metric_name}?per=day", token)
                update_daily_metric(metrics, metric_name, payload, item_key)
            except RuntimeError as exc:
                warnings.append(
                    f"Could not update {metric_name}: {exc}. "
                    "For traffic counters, configure TRAFFIC_TOKEN with repository Administration read access."
                )
    else:
        warnings.append("Could not update clone/view traffic: no TRAFFIC_TOKEN or GITHUB_TOKEN was available.")

    try:
        update_release_downloads(metrics, repository, token)
    except RuntimeError as exc:
        warnings.append(f"Could not update release asset downloads: {exc}.")

    metrics["last_updated"] = utc_now()
    metrics["warnings"] = warnings

    write_json(metrics_path, metrics)
    write_json(clones_badge_path, badge("recorded clones", int(metrics["clones"]["total"]), "blue"))
    write_json(
        releases_badge_path,
        badge("release downloads", int(metrics["release_asset_downloads"]["total"]), "brightgreen"),
    )

    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    print(f"Updated usage metrics for {repository}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
