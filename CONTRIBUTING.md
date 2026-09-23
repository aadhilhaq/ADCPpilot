# Contributing

Contributions should address a defined problem, explain the proposed change, and include appropriate verification. Describe changes to scientific behavior, defaults, or dependencies in the pull request.

## Reporting problems

Open a [GitHub issue](https://github.com/aadhilhaq/ADCP_docking/issues) or use **Send feedback** in the interface. Include:

- Software version or commit, operating system, and cluster environment.
- Steps to reproduce the problem and the expected result.
- Relevant settings and log excerpts, including the job's `.out` file when applicable.
- Output from `./check.sh` for installation problems.

Remove access tokens, credentials, and confidential data before sharing logs.

## Running tests

Run the full suite on an allocated compute node:

```bash
srun --pty -p cpu -c 8 -t 30 bash tests/run_all.sh
```

Adjust the partition and account for your cluster. The suite normally requires a Slurm allocation and reports each test file separately.

Tests are standalone scripts. Exit code `0` indicates a pass, `77` indicates a skip, and other codes indicate failure. A skipped test does not establish that the corresponding feature works. The suite does not replace a complete docking and rescoring run.

These individual tests do not require ADFRsuite or AmberTools:

```bash
python3 tests/test_api.py
python3 tests/test_engine.py
python3 tests/test_launcher.py
```

Include the test results and any skipped checks with your change.

## GitHub workflows

The test and package-publishing workflows are manually triggered. The desktop build workflow runs manually or when a tag matching `v*` is pushed. It builds Windows and macOS launchers on their respective platforms.

Tests that require external scientific programs or a cluster must also be checked in the appropriate environment. Keep reusable test logic in test scripts so it can be run locally as well as from a workflow.

## Adding tests

Follow the structure of the existing `tests/test_*.py` files. Report failures, use exit code `77` when required fixtures are unavailable, and print a concise result. Name each test for the behavior or failure it checks.

## Shared engine modules

`adcp_dock/_vendor/` contains modules copied from the HOPE peptide pipeline. Update these modules upstream and refresh the copies here to keep the implementations consistent.

```bash
python3 scripts/refresh_vendor.py --check
python3 scripts/refresh_vendor.py --audit
```

These commands check upstream differences and the vendor inventory. `tests/test_engine.py` checks integration with the shared modules.

## Code and documentation

Follow the existing Python structure and keep changes focused. Explain non-obvious decisions in comments. Update the documentation when commands, settings, dependencies, or outputs change.

Write documentation in direct language: state what the feature does, what the user must provide, and how to interpret the result. Distinguish implemented behavior from measured performance or scientific validation.

## Help

See [INSTALL.md](INSTALL.md) for installation and troubleshooting. For project questions, contact Aadhil Haq at aadhil.haq@tamu.edu or Dr. Sandun Fernando at sfernando@tamu.edu.
