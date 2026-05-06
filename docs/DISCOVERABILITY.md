# Discoverability And Packaging

This checklist is for maintainers who want ICLOCS to be easier to find, cite, install, and understand by humans, search engines, software indexes, and LMs.

## Repository Metadata

Recommended GitHub About metadata:

- Description: `MATLAB optimal-control software for direct collocation, trajectory optimization, and nonlinear programming workflows.`
- Website: `http://www.ee.ic.ac.uk/ICLOCS/`
- Topics: `matlab`, `optimal-control`, `direct-collocation`, `trajectory-optimization`, `nonlinear-programming`, `scientific-computing`, `ipopt`, `worhp`, `fmincon`, `imperial-college-london`.

Also add a social preview image in GitHub repository settings when a suitable project graphic is available.

## Machine-Readable Metadata

The repository includes:

- `CITATION.cff` for GitHub citation support.
- `codemeta.json` for research-software metadata indexes.
- `.zenodo.json` for Zenodo release metadata.
- `llms.txt` for LMs, agents, and repo-orientation tools.
- `docs/llms.txt` so the same orientation can be served from a future GitHub Pages docs site.

Review these files whenever the project version, citation, author list, website, or release process changes.

## Documentation Site

The `docs` folder includes a lightweight GitHub Pages landing page:

- `docs/index.md`
- `docs/_config.yml`

After this lands on `master`, maintainers can enable GitHub Pages from the repository settings and publish from the `docs` folder. If a future custom domain is used, update `README.md`, `codemeta.json`, `.zenodo.json`, and this page.

## Releases And Download Counts

Use GitHub Releases for user-visible release notes and durable downloads. Attach explicit assets when download counts matter:

- `ICLOCS-vX.Y.Z.zip`
- Optional `ICLOCS-vX.Y.Z.mltbx`

The automatic GitHub source-code archives are convenient, but explicit release assets are easier to name, describe, and count.

## External Indexes

Recommended external surfaces after a reviewed release:

- Zenodo GitHub integration for DOI-backed software releases.
- MATLAB File Exchange submission linked to the GitHub repository.
- The existing ICLOCS website should point to the GitHub repository, quick-install command, and release page.
- Any Imperial College software catalogue or research group page should use the same short description and keywords.

## Search Phrases To Cover

Keep these phrases in README/docs naturally, without keyword stuffing:

- MATLAB optimal control software.
- Direct collocation optimal control.
- Trajectory optimization in MATLAB.
- Nonlinear programming solver interface.
- IPOPT optimal control MATLAB.
- ICLOCS installation.
