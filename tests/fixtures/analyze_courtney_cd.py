#!/usr/bin/env python3
"""Compare measured free-flight Cd (Courtney et al. 2016, arXiv:1608.06500)
against the shapes of this engine's G1 and G7 standard drag curves.

WHAT THIS TESTS. A single ballistic coefficient asserts that a bullet's drag
curve is the standard curve scaled by a constant form factor:

    Cd_bullet(M) = i * Cd_standard(M),  i constant

That is the assumption behind quoting one BC for a bullet. This script fits the
best single `i` per bullet per model, then reports how far each measured point
still lands from that fit. The residual IS the error a constant-BC model makes,
isolated from muzzle velocity, sight height, twist and wind -- none of which
enter here.

A residual only means something if it exceeds the measurement's own SEM, which
is carried in the fixture and printed alongside.

Read-only: touches no engine code, only data/g1.csv and data/g7.csv.
"""
import csv, os, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SPEED_OF_SOUND_FPS = 1116.45   # 15 degC standard; see caveat in the summary


def load_curve(path):
    pts = []
    with open(path) as f:
        for row in csv.reader(f):
            if not row or row[0].strip().startswith("#"):
                continue
            try:
                pts.append((float(row[0]), float(row[1])))
            except ValueError:
                continue
    pts.sort()
    return pts


def interp(curve, mach):
    if mach <= curve[0][0]:
        return curve[0][1]
    if mach >= curve[-1][0]:
        return curve[-1][1]
    lo, hi = 0, len(curve) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if curve[mid][0] <= mach:
            lo = mid
        else:
            hi = mid
    (m0, c0), (m1, c1) = curve[lo], curve[hi]
    t = (mach - m0) / (m1 - m0)
    return c0 + t * (c1 - c0)


def load_fixture(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split(",")
            rows.append({
                "bullet": p[0], "dia": float(p[1]), "mass": float(p[2]),
                "v": float(p[3]), "cd": float(p[4]), "sem": float(p[5]),
                "n": int(p[6]), "table": p[8],
            })
    return rows


def main():
    g1 = load_curve(os.path.join(ROOT, "data", "g1.csv"))
    g7 = load_curve(os.path.join(ROOT, "data", "g7.csv"))
    rows = load_fixture(os.path.join(HERE, "courtney_2016_doppler_cd.csv"))
    print("  G1 curve: %d points, mach %.2f-%.2f" % (len(g1), g1[0][0], g1[-1][0]))
    print("  G7 curve: %d points, mach %.2f-%.2f" % (len(g7), g7[0][0], g7[-1][0]))

    groups = {}
    for r in rows:
        groups.setdefault((r["bullet"], r["dia"], r["mass"]), []).append(r)
    multi = {k: v for k, v in groups.items() if len(v) >= 2}

    print("\n  %d bullets have 2+ measured velocities (%d points total)\n"
          % (len(multi), sum(len(v) for v in multi.values())))

    summary = []
    for key in sorted(multi, key=lambda k: -len(multi[k])):
        bullet, dia, mass = key
        pts = sorted(multi[key], key=lambda r: -r["v"])
        machs = [p["v"] / SPEED_OF_SOUND_FPS for p in pts]
        print("=== %s %.3f\" %gGR   %d points, M%.2f-M%.2f ==="
              % (bullet, dia, mass, len(pts), min(machs), max(machs)))
        line = {}
        for name, curve in (("G1", g1), ("G7", g7)):
            std = [interp(curve, m) for m in machs]
            meas = [p["cd"] for p in pts]
            # best single form factor, least squares
            i = sum(s * c for s, c in zip(std, meas)) / sum(s * s for s in std)
            res = [(i * s - c) / c * 100.0 for s, c in zip(std, meas)]
            rms = math.sqrt(sum(r * r for r in res) / len(res))
            worst = max(res, key=abs)
            line[name] = (i, rms, worst)
            print("    %s  form factor i=%.3f   RMS residual %5.1f%%   worst %+.1f%%"
                  % (name, i, rms, worst))
        sem_typ = sum(p["sem"] for p in pts) / len(pts)
        print("    measurement SEM (mean of points): %.1f%%" % sem_typ)
        for p, m in zip(pts, machs):
            g1s, g7s = interp(g1, m), interp(g7, m)
            r1 = (line["G1"][0] * g1s - p["cd"]) / p["cd"] * 100
            r7 = (line["G7"][0] * g7s - p["cd"]) / p["cd"] * 100
            print("      M%.3f  measured Cd %.3f (SEM %.1f%%, n=%-2d)  G1 %+6.1f%%  G7 %+6.1f%%"
                  % (m, p["cd"], p["sem"], p["n"], r1, r7))
        summary.append((bullet, mass, dia, len(pts), sem_typ,
                        line["G1"][1], line["G7"][1]))
        print()

    print("=== summary: RMS residual of a CONSTANT form factor ===")
    print("  %-10s %6s %5s  %8s %8s %8s   %s"
          % ("bullet", "gr", "n_pts", "SEM%", "G1 RMS%", "G7 RMS%", "verdict"))
    for b, mass, dia, n, sem, r1, r7 in summary:
        best = "G7" if r7 < r1 else "G1"
        sig = "" if min(r1, r7) <= sem else "  <-- exceeds measurement error"
        print("  %-10s %6g %5d  %7.1f%% %7.1f%% %7.1f%%   %s%s"
              % (b, mass, n, sem, r1, r7, best, sig))


if __name__ == "__main__":
    main()
