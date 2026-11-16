#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""四方对照 h160 指标提取 + HTML 报告生成"""
import glob, json, os, re, sys

RUN = "/home/linfan/openfoam/run"
OUT_HTML = "/mnt/f/OpenFOAM/_cmp_h160_report.html"
NX, NY = 160, 320
PAPER = {"Vc1": 0.2472, "Vc2": 0.2315, "yc3": 1.1174}

CASES = [
    ("_cmp_h160_vofMULES", "VOF-MULES (incompressibleVoF)"),
    ("_cmp_h160_clsMULES", "CLSVOF-MULES (levelSetVoF)"),
    ("_cmp_h160_vofPLIC",  "VOF-PLIC (coupledAdvection/isoAdvector)"),
    ("_cmp_h160_clsPLIC",  "CLSVOF-PLIC (plicLevelSetVoF)"),
    ("_cmp_h160_clsPLIC_als", "CLSVOF-PLIC 强耦合 (advectNormal fromLevelSet)"),
]

def parse_dat(path):
    ts, vs = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("(")
            if len(parts) < 2:
                continue
            ts.append(float(parts[0].strip()))
            vs.append(float(parts[1].strip().rstrip(")").split()[1]))  # Uy
    return ts, vs

def parse_scalar_field(path):
    with open(path) as f:
        lines = f.read().splitlines()
    for i, ln in enumerate(lines):
        if not ln.startswith("internalField"):
            continue
        if ln.split()[1] != "nonuniform":
            return [float(ln.split()[-1])]
        j = i + 1
        try:
            n = int(ln.split()[-1])
        except ValueError:
            n = int(lines[j]); j += 1
        while lines[j].strip() != "(":
            j += 1
        j += 1
        vals = []
        while len(vals) < n:
            vals.append(float(lines[j])); j += 1
        return vals
    raise RuntimeError("internalField not found: " + path)

def yc_vol(path):
    alpha = parse_scalar_field(path)
    dx, dy, dz = 1.0/NX, 2.0/NY, 0.1
    num = den = 0.0
    for idx, a in enumerate(alpha):
        num += a * ((idx // NX) + 0.5) * dy
        den += a
    return (num/den if den > 0 else float("nan")), den*dx*dy*dz

def find_peaks(ts, vs, window=41, prominence=0.002, min_sep=0.3):
    n = len(ts); half = window//2
    sm = [sum(vs[max(0,i-half):min(n,i+half+1)])/(min(n,i+half+1)-max(0,i-half))
          for i in range(n)]
    peaks = []
    for i in range(1, n-1):
        if sm[i] > sm[i-1] and sm[i] > sm[i+1]:
            lo, hi = i, i
            while lo > 0 and sm[lo] <= sm[i] and sm[lo-1] <= sm[i]:
                lo -= 1
            while hi < n-1 and sm[hi] <= sm[i] and sm[hi+1] <= sm[i]:
                hi += 1
            if sm[i] - min(sm[lo], sm[hi]) < prominence:
                continue
            if peaks and ts[i]-peaks[-1][0] < min_sep:
                if sm[i] > peaks[-1][1]:
                    peaks[-1] = (ts[i], sm[i])
                continue
            peaks.append((ts[i], sm[i]))
    return peaks

def wall_time(case):
    p = os.path.join(RUN, case, "log.foamRun")
    try:
        m = re.findall(r"ClockTime = ([\d.]+) s", open(p).read())
        return float(m[-1]) if m else None
    except Exception:
        return None

results = []
for case, label in CASES:
    ts, vs = parse_dat(os.path.join(RUN, case, "postProcessing/bubbleMetrics/0/volFieldValue.dat"))
    peaks = find_peaks(ts, vs)
    ycs, vols = [], []
    dirs = sorted((d for d in glob.glob(os.path.join(RUN, case, "[0-9]*"))
                   if os.path.isfile(os.path.join(d, "alpha.air"))),
                  key=lambda x: float(os.path.basename(x)))
    for d in dirs:
        t = float(os.path.basename(d))
        yc, vol = yc_vol(os.path.join(d, "alpha.air"))
        ycs.append([t, yc]); vols.append([t, vol])
    vol0, volT = vols[0][1], vols[-1][1]
    drift = 100*(volT-vol0)/vol0 if vol0 else float("nan")
    res = {
        "case": case, "label": label,
        "ts": [round(t,6) for t in ts[::4]], "vc": [round(v,6) for v in vs[::4]],
        "yc_ts": [[round(t,4), round(y,5)] for t, y in ycs],
        "peaks": [[round(t,4), round(v,5)] for t, v in peaks[:4]],
        "Vc1": peaks[0][1] if len(peaks) > 0 else None,
        "Vc2": peaks[1][1] if len(peaks) > 1 else None,
        "yc3": ycs[-1][1],
        "vol0": vol0, "volT": volT, "drift": drift,
        "wall": wall_time(case),
    }
    results.append(res)
    print(f"{label}: Vc1={res['Vc1']:.4f} Vc2={res['Vc2']:.4f} "
          f"yc3={res['yc3']:.4f} drift={res['drift']:+.4f}% wall={res['wall']}s")

def dev(v, ref):
    return f"{100*(v-ref)/ref:+.2f}%"

# ---------------- HTML ----------------
colors = ["#e74c3c", "#27ae60", "#2980b9", "#8e44ad", "#e67e22"]
paper_row = ('<tr style="background:#f7f7f7"><td>Taran 论文 (1/h=160, VF-coupling)</td>'
             f'<td>{PAPER["Vc1"]:.4f}</td><td>{PAPER["Vc2"]:.4f}</td>'
             f'<td>{PAPER["yc3"]:.4f}</td><td>−0.0003%</td><td>—</td></tr>')
rows = ""
for r in results:
    rows += (f"<tr><td><b>{r['label']}</b></td>"
             f"<td>{r['Vc1']:.4f} ({dev(r['Vc1'],PAPER['Vc1'])})</td>"
             f"<td>{r['Vc2']:.4f} ({dev(r['Vc2'],PAPER['Vc2'])})</td>"
             f"<td>{r['yc3']:.4f} ({dev(r['yc3'],PAPER['yc3'])})</td>"
             f"<td>{r['drift']:+.4f}%</td>"
             f"<td>{r['wall']:.0f}</td></tr>")

chart_data = json.dumps({
    "series": [{"label": r["label"], "ts": r["ts"], "vc": r["vc"],
                "yc_ts": r["yc_ts"]} for r in results],
    "paper": PAPER,
})

html = f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>四方求解器对照 — Hysing 上升气泡 h160</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
<style>
body {{ font-family: "Segoe UI", "Microsoft YaHei", sans-serif; margin: 24px; background:#fff; color:#222; }}
h1 {{ font-size: 22px; }} h2 {{ font-size: 17px; margin-top: 28px; }}
table {{ border-collapse: collapse; margin-top: 8px; font-size: 14px; }}
th, td {{ border: 1px solid #ccc; padding: 6px 12px; text-align: center; }}
th {{ background: #2c3e50; color: #fff; }}
.chart {{ width: 100%; max-width: 900px; height: 380px; margin-top: 10px; }}
.note {{ color:#666; font-size: 13px; max-width: 900px; }}
.legend {{ margin-top: 4px; font-size: 13px; color:#444; }}
</style></head><body>
<h1>四方求解器对照 — Hysing 上升气泡基准（h160，160×320）</h1>
<p class="note">同一网格 / 同一初始场（逐字节一致）/ 同一固定 dt=3.90625e-04 / endTime 3s，唯一变量为求解器。
Vc,1、Vc,2 为气泡平均上升速度曲线的前两个显著峰值（t≈1s、t≈2s），yc(3s) 为 3s 时刻气泡质心高度，漂移为气泡体积变化。
参考行：Taran et al. 论文 Table 2（1/h=160，volume-fraction coupling）。</p>

<h2>指标对照表</h2>
<table><tr><th>求解器</th><th>Vc,1</th><th>Vc,2</th><th>yc(3s)</th><th>体积漂移</th><th>墙钟 (s)</th></tr>
{paper_row}
{rows}
</table>

<h2>Vc(t) — 气泡平均上升速度</h2>
<div class="chart"><canvas id="c1"></canvas></div>
<p class="legend">虚线为 Taran 论文参考值不可直接叠加曲线，仅作对照：Vc,1={PAPER['Vc1']}, Vc,2={PAPER['Vc2']}。</p>

<h2>yc(t) — 气泡质心高度</h2>
<div class="chart"><canvas id="c2"></canvas></div>

<script>
const D = {chart_data};
const colors = {json.dumps(colors)};
new Chart(document.getElementById('c1'), {{
  type: 'line',
  data: {{ datasets: D.series.map((s,i)=>({{
    label: s.label, data: s.ts.map((t,k)=>({{x:t, y:s.vc[k]}})),
    borderColor: colors[i], backgroundColor: colors[i],
    borderWidth: 1.6, pointRadius: 0, tension: 0.15, parsing: false
  }})) }},
  options: {{
    parsing: false,
    scales: {{ x: {{ type:'linear', title:{{display:true,text:'t (s)'}} }},
              y: {{ title:{{display:true,text:'Vc (m/s)'}} }} }},
    plugins: {{ legend: {{ position:'top' }}, interaction: {{ mode:'nearest', intersect:false }} }},
    elements: {{ point: {{ radius: 0 }} }}
  }}
}});
new Chart(document.getElementById('c2'), {{
  type: 'line',
  data: {{ datasets: D.series.map((s,i)=>({{
    label: s.label, data: s.yc_ts.map(p=>({{x:p[0], y:p[1]}})),
    borderColor: colors[i], borderWidth: 1.6, pointRadius: 0, tension: 0.15
  }})) }},
  options: {{
    parsing: false,
    scales: {{ x: {{ type:'linear', title:{{display:true,text:'t (s)'}} }},
              y: {{ title:{{display:true,text:'yc (m)'}} }} }},
    plugins: {{ legend: {{ position:'top' }} }},
    elements: {{ point: {{ radius: 0 }} }}
  }}
}});
</script>
</body></html>"""

with open(OUT_HTML, "w", encoding="utf-8") as f:
    f.write(html)
print("HTML ->", OUT_HTML)
