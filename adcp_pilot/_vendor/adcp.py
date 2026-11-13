"""ADCP's ranked poses, read for everything downstream.

By Aadhil Haq

ADCP writes each ranked pose with all its atoms, side chains included, as
out_ranked_<n>.pdb (split from its multi-model out_out.pdb by
adcp11.split_poses). Its replicas' own trajectories are in ADCP's reduced
representation, with side chains cut at CG and no aromatic ring atoms, so
contacts counted on them are wrong; nothing here reads them.
"""
from __future__ import annotations
import os
import numpy as np
from .receptor import is_hydrogen


def read_pose(wd: str, rank: int = 1) -> dict | None:
    """The all-atom ranked pose, without hydrogens, or None."""
    p = os.path.join(wd, f"out_ranked_{rank}.pdb")
    if not os.path.isfile(p):
        return None
    xyz, resnum, resname, name = [], [], [], []
    for l in open(p):
        if not l.startswith("ATOM"):
            continue
        nm = l[12:16].strip()
        el = l[76:78].strip()
        if is_hydrogen(nm, el):
            continue
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        resnum.append(int(l[22:26])); resname.append(l[17:20].strip()); name.append(nm)
    if not xyz:
        return None
    return dict(xyz=np.array(xyz), resnum=np.array(resnum),
                resname=resname, name=name)
