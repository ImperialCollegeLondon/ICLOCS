# ICLOCS
Imperial College London Optimal Control Software (ICLOCS)

[![Recorded clones](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/ImperialCollegeLondon/ICLOCS/usage-metrics/docs/metrics/clones-badge.json)](https://github.com/ImperialCollegeLondon/ICLOCS/blob/usage-metrics/docs/metrics/github-usage.json)
[![Release downloads](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/ImperialCollegeLondon/ICLOCS/usage-metrics/docs/metrics/release-downloads-badge.json)](https://github.com/ImperialCollegeLondon/ICLOCS/blob/usage-metrics/docs/metrics/github-usage.json)
[![MATLAB CI](https://github.com/ImperialCollegeLondon/ICLOCS/actions/workflows/matlab-ci.yml/badge.svg)](https://github.com/ImperialCollegeLondon/ICLOCS/actions/workflows/matlab-ci.yml)

http://www.ee.ic.ac.uk/ICLOCS/

Modernization work is tracked in [`TODO.md`](TODO.md).

## Documentation

- [Installation](docs/INSTALLATION.md)
- [Examples](docs/EXAMPLES.md)
- [Example results](docs/EXAMPLE_RESULTS.md)
- [Problem formulation](docs/PROBLEM_FORMULATION.md)
- [Options reference](docs/OPTIONS_REFERENCE.md)
- [Closed-loop simulation](docs/CLOSED_LOOP_SIMULATION.md)
- [Development and CI](docs/DEVELOPMENT.md)
- [API index](docs/API-INDEX.md)
- [Repository layout](docs/REPOSITORY_LAYOUT.md)
- [Release guide](docs/RELEASE_GUIDE.md)
- [Source website map](docs/SOURCE_WEBSITE.md)

## Repository usage counter

This repository includes a lightweight GitHub Actions workflow that records repository usage metrics in `docs/metrics/github-usage.json` on the `usage-metrics` branch. The README badges above are updated from the same metric files once the workflow has run on GitHub.

GitHub exposes clone traffic only for the previous 14 days, so the recorded clone counter starts from the first successful scheduled or manual workflow run and cannot reconstruct historical downloads. GitHub release asset downloads are counted separately when releases provide downloadable assets. To enable clone/view traffic collection, add a repository secret named `TRAFFIC_TOKEN` with repository Administration read access, then run the **Update usage metrics** workflow from the Actions tab.

## Note for updating your problem from version 2 to run with version 2.5

We have recently introduced changes to the problem file (see myProblem.m in the template). 

To update your existing problem to work with version 2.5.1 onwards, please make the following changes:
* Update the codes for section "Get function handles and return to Main.m" with the ones in the new template
* Remove all codes below the divider

	%% Leave the following unchanged! %%
