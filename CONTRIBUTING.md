# Contributing

## Reporting problems

Open a [GitHub issue](https://github.com/aadhilhaq/ADCP_docking/issues) or use **Send feedback** in the interface. Include the version or commit, the operating system and cluster, the steps that reproduce the problem, the settings, and the relevant log (the job's `.out` file where there is one). For an installation problem, include the output of `adcp-pilot check`. Remove access tokens, passwords and confidential data first.

## Running the tests

On a compute node, adding the account and partition your cluster requires:

```bash
srun --pty -c 8 -t 30 bash tests/run_all.sh
```

Each test file is a standalone script: exit code 0 is a pass, 77 a skip (a program it needs is missing), anything else a failure. A skipped test has checked nothing. These need neither the ADCP 1.1 suite nor AmberTools:

```bash
python3 tests/test_api.py
python3 tests/test_engine.py
python3 tests/test_launcher.py
```

## Changes

Keep a change focused. Say what it changes in the scientific behavior, the defaults or the dependencies, update the documentation when a command, setting or output changes, and include the test results.

## Contact

Aadhil Haq (aadhil.haq@tamu.edu) or Dr. Sandun Fernando (sfernando@tamu.edu).
