#!/usr/bin/env python3
"""Exact cut-cell alpha.air for a static droplet benchmark.

Circle radius r centred at (xc, yc) in a W x H x 0.1 domain, uniform
Cartesian mesh nx x ny x 1. alpha.air = chord-area fraction via
Gauss-Legendre quadrature with sub-interval splitting at tangency points.

Usage: drop_alpha.py <nx> <ny> <W> <H> <xc> <yc> <r> > 0/alpha.air
"""
import sys
import math

NQ = 24  # Gauss-Legendre points per smooth sub-interval


def gauss_legendre(n):
    nodes, weights = [], []
    # numpy-free composite: use midpoint fallback if numpy missing
    try:
        import numpy as np
        x, w = np.polynomial.legendre.leggauss(n)
        return list(x), list(w)
    except ImportError:
        step = 2.0 / n
        for i in range(n):
            nodes.append(-1 + (i + 0.5) * step)
            weights.append(step)
        return nodes, weights


def chord_width(y, xc, yc, r, x0, x1):
    dy = y - yc
    s2 = r * r - dy * dy
    if s2 <= 0.0:
        return 0.0
    s = math.sqrt(s2)
    lo = max(x0, xc - s)
    hi = min(x1, xc + s)
    return hi - lo if hi > lo else 0.0


def cell_alpha(x0, x1, y0, y1, xc, yc, r):
    # split y-interval at tangency heights so integrand is smooth
    cuts = {y0, y1}
    if abs(yc - yc) >= 0:
        for s in (yc - r, yc + r):
            if y0 < s < y1:
                cuts.add(s)
    cuts = sorted(cuts)
    total = 0.0
    xs, ws = gauss_legendre(NQ)
    for a, b in zip(cuts[:-1], cuts[1:]):
        for x, w in zip(xs, ws):
            y = 0.5 * (a + b) + 0.5 * (b - a) * x
            total += w * chord_width(y, xc, yc, r, x0, x1)
    total *= 0.5 * (b - a) if False else 0.0
    # redo properly: accumulate per sub-interval
    total = 0.0
    for a, b in zip(cuts[:-1], cuts[1:]):
        sub = 0.0
        for x, w in zip(xs, ws):
            y = 0.5 * (a + b) + 0.5 * (b - a) * x
            sub += w * chord_width(y, xc, yc, r, x0, x1)
        total += 0.5 * (b - a) * sub
    area = total
    return area / ((x1 - x0) * (y1 - y0))


def main():
    nx, ny = int(sys.argv[1]), int(sys.argv[2])
    W, H = float(sys.argv[3]), float(sys.argv[4])
    xc, yc, r = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
    dx, dy = W / nx, H / ny
    print("FoamFile")
    print("{")
    print("    format      ascii;")
    print("    class       volScalarField;")
    print("    location    \"0\";")
    print("    object      alpha.air;")
    print("}")
    print("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")
    print("")
    print("dimensions      [0 0 0 0 0 0 0];")
    print("")
    print("internalField   nonuniform List<scalar>")
    print(f"{nx * ny}")
    print("(")
    for j in range(ny):
        y0, y1 = j * dy, (j + 1) * dy
        for i in range(nx):
            x0, x1 = i * dx, (i + 1) * dx
            a = cell_alpha(x0, x1, y0, y1, xc, yc, r)
            print(f"    {a:.12f}")
    print(")")
    print(";")
    print("")
    print("boundaryField")
    print("{")
    print("    bottom { type zeroGradient; }")
    print("    top    { type zeroGradient; }")
    print("    walls  { type zeroGradient; }")
    print("    frontAndBack { type empty; }")
    print("}")
    # self-check
    exact = math.pi * r * r
    vol = sum(
        cell_alpha(i * dx, (i + 1) * dx, j * dy, (j + 1) * dy, xc, yc, r) * dx * dy
        for j in range(ny)
        for i in range(nx)
    )
    print(f"// volume check: {vol:.12f} vs exact {exact:.12f} "
          f"rel.err {abs(vol - exact) / exact:.2e}", file=sys.stderr)


if __name__ == "__main__":
    main()
