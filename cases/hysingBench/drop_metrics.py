#!/usr/bin/env python3
"""Spurious-currents metrics for the static droplet benchmark.

For each written time of a case:
  - max mag(U) over all cells (parasite velocity magnitude)
  - band statistics of curvature-like fields if written
      (levelSetCurvature for LS solvers, K for interfaceProperties if written)
  - droplet volume drift from alpha.air
Usage: drop_metrics.py <case>
"""
import glob
import math
import os
import re
import sys


def parse_vector_field(path):
    """Return list of [x,y,z] internalField values of an ascii vector field."""
    with open(path) as f:
        lines = f.read().splitlines()
    vals = []
    for i, ln in enumerate(lines):
        if ln.startswith("internalField"):
            toks = ln.split()
            if len(toks) > 1 and toks[1] == "uniform":
                m = re.findall(r"[-+0-9.eE]+", ln)
                vals = [[float(m[-3]), float(m[-2]), float(m[-1])]]
                return vals
            # nonuniform: find count then read "( ... )"
            n = None
            rest = ln.split("nonuniform")
            m = re.search(r"nonuniform\s+List<\w+>\s*(\d+)", ln)
            if m:
                n = int(m.group(1))
            else:
                j = i + 1
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j < len(lines):
                    try:
                        n = int(lines[j].strip())
                    except ValueError:
                        return vals
            j = i + 1
            # find opening paren
            while j < len(lines) and "(" not in lines[j]:
                j += 1
            buf = []
            while j < len(lines) and len(buf) < (n or 0):
                for mm in re.finditer(r"\(([^)]*)\)", lines[j]):
                    if len(buf) < (n or 0):
                        comp = mm.group(1).split()
                        if len(comp) == 3:
                            buf.append([float(c) for c in comp])
                j += 1
            vals = buf
            return vals
    return vals


def parse_scalar_field(path):
    with open(path) as f:
        lines = f.read().splitlines()
    vals = []
    for i, ln in enumerate(lines):
        if ln.startswith("internalField"):
            m = re.search(r"nonuniform\s+List<scalar>\s*(\d+)", ln)
            n = None
            if m:
                n = int(m.group(1))
            else:
                # v14 style: count may sit on the next line
                j = i + 1
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j < len(lines) and lines[j].strip().isdigit():
                    n = int(lines[j].strip())
                else:
                    m2 = re.search(r"uniform\s+\(?([-\d.eE+]+)\)?", ln)
                    if m2:
                        return [float(m2.group(1))]
                    return vals
            j = i + 1
            while j < len(lines) and "(" not in lines[j]:
                j += 1
            j += 1
            buf = []
            while j < len(lines) and len(buf) < n:
                for mm in re.finditer(r"[-+0-9.eE]+", lines[j]):
                    if len(buf) < n:
                        buf.append(float(mm.group(0)))
                j += 1
            vals = buf
            return vals
    return vals


def main():
    case = sys.argv[1]
    print(f"===== {case} =====")
    print(f"{'t':>10s}  {'max|U|':>12s}  {'Vdrop':>12s}")
    tdirs = sorted(
        (d for d in glob.glob(os.path.join(case, "[0-9]*")) if os.path.isdir(d)),
        key=lambda d: float(os.path.basename(d)),
    )
    v0 = None
    for td in tdirs:
        ufile = os.path.join(td, "U")
        afile = os.path.join(td, "alpha.air")
        t = float(os.path.basename(td))
        if not os.path.isfile(ufile):
            continue
        us = parse_vector_field(ufile)
        mags = [math.sqrt(x * x + y * y + z * z) for x, y, z in us]
        mx = max(mags) if mags else float("nan")
        vol = None
        if os.path.isfile(afile):
            al = parse_scalar_field(afile)
            n = len(al)
            dx = dy = 1.0 / 80.0
            vol = sum(a for a in al) * dx * dy * 0.1
            if v0 is None:
                v0 = vol
        drift = (vol - v0) / v0 * 100 if (vol is not None and v0) else float("nan")
        print(f"{t:10.3f}  {mx:12.4e}  {drift:+11.4f}%")


if __name__ == "__main__":
    main()
