"""A finished run to open in the interface before docking anything.

MDM2 (PDB 1YCR, chain A) with three peptides: p53 residues 17 to 29, PMI,
and RFMDYWEGL, the core of the phage display peptide 12/1. Docked with ADCP
(25 replicas of 2,000,000 steps) and rescored with MM-GBSA on five poses.

The run is stored in the package with its own folder written as @RUN@, and
copied into a runs folder on request:

    adcp-dock example [--runs DIR]
"""
import json
import os
import shutil

NAME = "example_MDM2"
PLACEHOLDER = "@RUN@"
SOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "example_run")


def _fill(value, dest):
    """Every string in a JSON value with the placeholder replaced by dest."""
    if isinstance(value, str):
        return value.replace(PLACEHOLDER, dest)
    if isinstance(value, list):
        return [_fill(v, dest) for v in value]
    if isinstance(value, dict):
        return {k: _fill(v, dest) for k, v in value.items()}
    return value


def install(runs, force=False):
    """Copy the example into the runs folder and return where it went.

    An example that is already there is left as it is, unless force is set,
    so opening the interface again keeps whatever was done to it.
    """
    runs = os.path.abspath(os.path.expanduser(runs))
    dest = os.path.join(runs, NAME)
    if os.path.exists(dest):
        if not force:
            return dest
        shutil.rmtree(dest)
    os.makedirs(runs, exist_ok=True)
    shutil.copytree(SOURCE, dest)
    for base, _dirs, names in os.walk(dest):
        for n in names:
            path = os.path.join(base, n)
            if n.endswith(".json"):
                # Through the JSON parser rather than as text, so that a
                # folder name with a quote or a backslash stays valid JSON.
                with open(path, encoding="utf-8") as fh:
                    data = json.load(fh)
                with open(path, "w", encoding="utf-8") as fh:
                    json.dump(_fill(data, dest), fh, indent=1)
            elif n.endswith((".txt", ".tsv", ".out", ".log", ".bib", ".ris")):
                with open(path, encoding="utf-8", errors="replace") as fh:
                    text = fh.read()
                if PLACEHOLDER in text:
                    with open(path, "w", encoding="utf-8") as fh:
                        fh.write(text.replace(PLACEHOLDER, dest))
    return dest
