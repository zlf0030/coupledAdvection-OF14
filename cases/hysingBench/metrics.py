#!/usr/bin/env python3
"""Extract Hysing case-2 benchmark quantities from a completed baseline run.

Vc(t): from postProcessing/bubbleMetrics/0/volFieldValue.dat (per-step
       volAverage of U weighted by alpha.air). Vc = Uy component.
yc(t): from written alpha.air fields + analytic uniform-mesh cell centres.

Reports Vc,1 / Vc,2 (first two velocity peaks and their times) and yc(t=3),
compared with the paper's volume-fraction-coupling row (1/h=160):
    Vc,1 = 0.2472  Vc,2 = 0.2315  yc(3s) = 1.1174
"""
import glob
import math
import os
import sys

# ---------------------------------------------------------------- reference
PAPER = {"Vc1": 0.2472, "Vc2": 0.2315, "yc3": 1.1174}


def parse_dat(path):
    ts, vs = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("(")
            if len(parts) < 2:
                continue
            t = float(parts[0].strip())
            comp = parts[1].strip().rstrip(")").split()
            ts.append(t)
            vs.append(float(comp[1]))  # Uy
    return ts, vs


def parse_scalar_field(path):
    """Parse internalField (nonuniform or uniform) of an ascii field file."""
    with open(path) as f:
        lines = f.read().splitlines()
    for i, ln in enumerate(lines):
        if not ln.startswith("internalField"):
            continue
        if ln.split()[1] != "nonuniform":
            return [float(ln.split()[-1])]
        # nonuniform: count may be on this line or the next
        rest = ln.split()[1:]
        j = i + 1
        try:
            n = int(rest[-1])
        except ValueError:
            n = int(lines[j])
            j += 1
        while lines[j].strip() != "(":
            j += 1
        j += 1
        vals = []
        while len(vals) < n:
            vals.append(float(lines[j]))
            j += 1
        return vals
    raise RuntimeError(f"internalField not found in {path}")


def yc_from_field(path, nx, ny):
    alpha = parse_scalar_field(path)
    dx, dy = 1.0 / nx, 2.0 / ny
    dz = 0.1
    num = 0.0
    den = 0.0
    for idx, a in enumerate(alpha):
        j = idx // nx
        ycen = (j + 0.5) * dy
        num += a * ycen
        den += a
    return num / den if den > 0 else float("nan"), den * dx * dy * dz


def find_peaks(ts, vs, window=41, prominence=0.002, min_sep=0.3):
    n = len(ts)
    half = window // 2
    sm = []
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        sm.append(sum(vs[lo:hi]) / (hi - lo))
    peaks = []
    for i in range(1, n - 1):
        if sm[i] > sm[i - 1] and sm[i] > sm[i + 1]:
            # prominence vs neighbouring valleys
            lo, hi = i, i
            while lo > 0 and sm[lo] <= sm[i] and sm[lo - 1] <= sm[i]:
                lo -= 1
            while hi < n - 1 and sm[hi] <= sm[i] and sm[hi + 1] <= sm[i]:
                hi += 1
            base = min(sm[lo], sm[hi])
            if sm[i] - base >= prominence:
                if peaks and ts[i] - peaks[-1][0] < min_sep:
                    if sm[i] > peaks[-1][1]:
                        peaks[-1] = (ts[i], sm[i])
                    continue
                peaks.append((ts[i], sm[i]))
    return peaks


def main():
    case = sys.argv[1]
    nx, ny = int(sys.argv[2]), int(sys.argv[3])
    print(f"===== {case} (1/h = {1.0/(1.0/nx):.0f}) =====")

    dat = os.path.join(case, "postProcessing/bubbleMetrics/0/volFieldValue.dat")
    ts, vs = parse_dat(dat)
    print(f"Vc(t): {len(ts)} 个时间点, {ts[0]:.4f} -> {ts[-1]:.4f} s")

    peaks = find_peaks(ts, vs)
    print("Vc(t) 局部极大值（前 4 个）:")
    for t, v in peaks[:4]:
        print(f"    t = {t:.4f}  Vc = {v:.4f}")

    # yc at written times
    ycs = []
    for d in sorted(glob.glob(os.path.join(case, "[0-9]*")),
                    key=lambda x: float(os.path.basename(x))):
        fp = os.path.join(d, "alpha.air")
        if not os.path.isfile(fp):
            continue
        t = float(os.path.basename(d))
        try:
            yc, vol = yc_from_field(fp, nx, ny)
            ycs.append((t, yc, vol))
        except Exception as e:
            print(f"  (skip t={t}: {e})")
    if ycs:
        vol0, volT = ycs[0][2], ycs[-1][2]
        print(f"yc(t): {len(ycs)} 个时刻")
        # 打印若干样本
        for t, yc, _ in ycs[::max(1, len(ycs)//8)][1:]:
            print(f"    t = {t:.3f}  yc = {yc:.4f}")
        print(f"    t = {ycs[-1][0]:.3f}  yc = {ycs[-1][1]:.4f}   <-- yc(t_end)")
        if vol0:
            print(f"气泡体积: t=0 {vol0:.8f}  t_end {volT:.8f}  "
                  f"变化 {100*(volT-vol0)/vol0:+.4f}%")

    print("\n--- 与论文 volume-fraction coupling 行 (1/h=160) 对比 ---")
    rows = []
    if len(peaks) >= 1:
        rows.append(("Vc,1", peaks[0][1], PAPER["Vc1"]))
    if len(peaks) >= 2:
        rows.append(("Vc,2", peaks[1][1], PAPER["Vc2"]))
    if ycs:
        rows.append(("yc(3s)", ycs[-1][1], PAPER["yc3"]))
    for name, val, ref in rows:
        dev = 100 * (val - ref) / ref
        flag = "PASS" if abs(dev) <= 5 else "FAIL"
        print(f"    {name:8s} local = {val:.4f}  paper = {ref:.4f}  "
              f"偏差 {dev:+.2f}%  [{flag}]")


if __name__ == "__main__":
    main()
