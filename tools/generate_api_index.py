"""Generate a Markdown inventory of MATLAB files and first help lines."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs" / "API-INDEX.md"
MATLAB_ROOTS = [
    "src",
    "usr",
    "exampleProblems",
    "exampleProblems_fmincon",
    "closedLoopExamples",
    "tools",
]


def first_help_line(path: Path) -> tuple[str, str]:
    kind = "script"
    help_line = ""

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("function"):
            kind = "function"
            continue
        if line.startswith("%"):
            text = line.lstrip("%").strip()
            if text and set(text) != {"-"}:
                help_line = text
                break
            continue
        if help_line:
            break
        if not line.startswith("%"):
            break

    return kind, help_line or "Needs a MATLAB help summary."


def iter_matlab_files() -> list[Path]:
    files: list[Path] = []
    for root in MATLAB_ROOTS:
        directory = ROOT / root
        if directory.exists():
            files.extend(directory.rglob("*.m"))
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix().lower())


def main() -> int:
    rows = []
    missing = 0
    for path in iter_matlab_files():
        kind, summary = first_help_line(path)
        if summary == "Needs a MATLAB help summary.":
            missing += 1
        rel_path = path.relative_to(ROOT).as_posix()
        rows.append(f"| `{rel_path}` | {kind} | {summary} |")

    content = [
        "# API Index",
        "",
        "This generated inventory lists MATLAB files and their first help summary line. It is a documentation triage aid, not a replacement for full function-level documentation.",
        "",
        f"- MATLAB files indexed: {len(rows)}",
        f"- Files needing a first help summary: {missing}",
        "",
        "Regenerate with:",
        "",
        "```bash",
        "python3 tools/generate_api_index.py",
        "```",
        "",
        "| File | Type | Current help summary |",
        "| --- | --- | --- |",
        *rows,
        "",
    ]
    OUTPUT.write_text("\n".join(content), encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
