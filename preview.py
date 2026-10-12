#!/usr/bin/env python3
"""Draw SVG previews (cross-section and level plan) from a generated dataset.

    python3 preview.py --data data --out docs
"""
import argparse
import csv
import os
from collections import defaultdict

from porphyry_gen.model import E0, N0

LITHO_FILL = {"QMP": "#e8a3b5", "IMP": "#f2c48d", "BXH": "#7b5aa6", "DAC": "#8fb09a"}
ALT_FILL = {"PRO": "#9fd39a", "POT": "#f0a080", "PHY": "#f3e27a", "ARG": "#b9d7f0"}
SHELL_STROKE = {"CU_0.2": "#2a6fdb", "CU_0.4": "#d98a00", "CU_0.8": "#c0282d"}
LINE_STROKE = {"TOPO": "#222", "BOC": "#8a6d3b", "BOLC": "#b5651d", "BOX": "#2e8b57", "TOH": "#1f4e9c"}
CU_BINS = [(0.1, "#c7c7c7"), (0.2, "#5aa0e6"), (0.4, "#3bb273"), (0.8, "#f0a020"), (99, "#d62728")]


def read_strings(path):
    out = defaultdict(list)  # (view_label, unit) -> list of (closed, pts)
    cur, key, closed = None, None, 0
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["STRING_ID"] != cur:
                cur = r["STRING_ID"]
                key = (r["SECTION_OR_PLAN"], r["UNIT"])
                closed = int(r["CLOSED"])
                out[key].append((closed, []))
            out[key][-1][1].append((float(r["EAST"]), float(r["NORTH"]), float(r["RL"])))
    return out


def cu_colour(v):
    for lim, col in CU_BINS:
        if v < lim:
            return col
    return CU_BINS[-1][1]


class Panel:
    def __init__(self, x0, y0, w, h, umin, umax, vmin, vmax):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.umin, self.umax, self.vmin, self.vmax = umin, umax, vmin, vmax
        self.parts = []

    def xy(self, u, v):
        return (self.x0 + (u - self.umin) / (self.umax - self.umin) * self.w,
                self.y0 + (self.vmax - v) / (self.vmax - self.vmin) * self.h)

    def path(self, pts, closed, **style):
        d = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(self.xy(u, v) for u, v in pts))
        st = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in style.items())
        self.parts.append(f'<path d="{d}{" Z" if closed else ""}" {st}/>')

    def frame(self, title, ulabel, vlabel, step):
        p = [f'<rect x="{self.x0}" y="{self.y0}" width="{self.w}" height="{self.h}" fill="none" stroke="#333"/>',
             f'<text x="{self.x0}" y="{self.y0 - 10}" font-size="15" font-weight="bold">{title}</text>']
        u = (self.umin // step + 1) * step
        while u < self.umax:
            x, _ = self.xy(u, self.vmin)
            p.append(f'<text x="{x:.0f}" y="{self.y0 + self.h + 14}" font-size="10" text-anchor="middle">{u:.0f}</text>')
            u += step
        v = (self.vmin // step + 1) * step
        while v < self.vmax:
            _, y = self.xy(self.umin, v)
            p.append(f'<text x="{self.x0 - 4}" y="{y + 3:.0f}" font-size="10" text-anchor="end">{v:.0f}</text>')
            v += step
        p.append(f'<text x="{self.x0 + self.w / 2}" y="{self.y0 + self.h + 30}" font-size="11" text-anchor="middle">{ulabel}</text>')
        p.append(f'<text x="{self.x0 - 48}" y="{self.y0 + self.h / 2}" font-size="11" text-anchor="middle" '
                 f'transform="rotate(-90 {self.x0 - 48} {self.y0 + self.h / 2})">{vlabel}</text>')
        clip = f"clip{id(self)}"
        return (f'<clipPath id="{clip}"><rect x="{self.x0}" y="{self.y0}" width="{self.w}" height="{self.h}"/></clipPath>'
                f'<g clip-path="url(#{clip})">' + "".join(self.parts) + "</g>" + "".join(p))


def legend(x, y, items, kind="rect"):
    out = []
    for i, (label, col) in enumerate(items):
        yy = y + i * 17
        if kind == "rect":
            out.append(f'<rect x="{x}" y="{yy}" width="14" height="11" fill="{col}" stroke="#555" stroke-width="0.5"/>')
        elif kind == "line":
            out.append(f'<line x1="{x}" y1="{yy + 6}" x2="{x + 14}" y2="{yy + 6}" stroke="{col}" stroke-width="2"/>')
        else:
            out.append(f'<circle cx="{x + 7}" cy="{yy + 6}" r="4" fill="{col}"/>')
        out.append(f'<text x="{x + 20}" y="{yy + 10}" font-size="11">{label}</text>')
    return "".join(out)


def section_svg(data, strings, label, north):
    samples = []
    with open(os.path.join(data, "drillholes", "samples_desurveyed.csv")) as fh:
        for r in csv.DictReader(fh):
            if abs(float(r["NORTH"]) - north) < 50 and r["CuT_pct"]:
                samples.append((float(r["EAST"]), float(r["RL"]), float(r["CuT_pct"])))
    W, H = 1180, 1060
    ext = (E0 - 1000, E0 + 1000, 300, 1420)
    a = Panel(70, 50, 860, 440, *ext)
    b = Panel(70, 570, 860, 440, *ext)
    get = lambda unit: strings.get((label, unit), [])
    for unit, col in LITHO_FILL.items():
        for closed, pts in get(unit):
            a.path([(p[0], p[2]) for p in pts], closed, fill=col, stroke="#444", stroke_width=0.6)
    for unit, col in LINE_STROKE.items():
        for _, pts in get(unit):
            a.path([(p[0], p[2]) for p in pts], False, fill="none", stroke=col, stroke_width=1.6 if unit == "TOPO" else 1.1,
                   stroke_dasharray="" if unit == "TOPO" else "6 3")
    for e, z, cu in samples:
        x, y = a.xy(e, z)
        a.parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.6" fill="{cu_colour(cu)}"/>')
    for unit, col in ALT_FILL.items():
        for closed, pts in get(unit):
            b.path([(p[0], p[2]) for p in pts], closed, fill=col, fill_opacity=0.85, stroke="#666", stroke_width=0.5)
    for unit, col in LITHO_FILL.items():
        for closed, pts in get(unit):
            b.path([(p[0], p[2]) for p in pts], closed, fill="none", stroke="#555", stroke_width=0.6, stroke_dasharray="2 2")
    for unit, col in SHELL_STROKE.items():
        for closed, pts in get(unit):
            b.path([(p[0], p[2]) for p in pts], closed, fill="none", stroke=col, stroke_width=1.8)
    for _, pts in get("TOPO"):
        b.path([(p[0], p[2]) for p in pts], False, fill="none", stroke="#222", stroke_width=1.6)
    body = a.frame(f"Section {label} — geology, weathering profile and drill samples (±50 m)", "Easting (m)", "RL (m)", 200)
    body += b.frame(f"Section {label} — alteration and CuT grade-shell interpretation", "Easting (m)", "RL (m)", 200)
    lx = 950
    body += '<text x="950" y="62" font-size="12" font-weight="bold">Lithology</text>'
    body += legend(lx, 70, [("QMP stock", LITHO_FILL["QMP"]), ("IMP dyke", LITHO_FILL["IMP"]), ("BXH breccia", LITHO_FILL["BXH"]),
                            ("DAC dyke", LITHO_FILL["DAC"]), ("AND host (blank)", "#ffffff")])
    body += '<text x="950" y="172" font-size="12" font-weight="bold">Surfaces</text>'
    body += legend(lx, 180, [("Topography", LINE_STROKE["TOPO"]), ("Base of cover", LINE_STROKE["BOC"]),
                             ("Base leached cap", LINE_STROKE["BOLC"]), ("Base of oxide", LINE_STROKE["BOX"]),
                             ("Top of hypogene", LINE_STROKE["TOH"])], "line")
    body += '<text x="950" y="282" font-size="12" font-weight="bold">CuT % samples</text>'
    body += legend(lx, 290, [("&lt; 0.1", CU_BINS[0][1]), ("0.1 – 0.2", CU_BINS[1][1]), ("0.2 – 0.4", CU_BINS[2][1]),
                             ("0.4 – 0.8", CU_BINS[3][1]), ("≥ 0.8", CU_BINS[4][1])], "dot")
    body += '<text x="950" y="582" font-size="12" font-weight="bold">Alteration</text>'
    body += legend(lx, 590, [("Argillic", ALT_FILL["ARG"]), ("Phyllic (QSP)", ALT_FILL["PHY"]), ("Potassic", ALT_FILL["POT"]),
                             ("Propylitic", ALT_FILL["PRO"])])
    body += '<text x="950" y="675" font-size="12" font-weight="bold">CuT grade shells</text>'
    body += legend(lx, 683, [("0.2 %", SHELL_STROKE["CU_0.2"]), ("0.4 %", SHELL_STROKE["CU_0.4"]), ("0.8 %", SHELL_STROKE["CU_0.8"])], "line")
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
            f'font-family="Helvetica, Arial, sans-serif"><rect width="100%" height="100%" fill="#fff"/>{body}</svg>')


def plan_svg(data, strings, label):
    W, H = 1000, 860
    p = Panel(70, 50, 780, 780, E0 - 1000, E0 + 1000, N0 - 1000, N0 + 1000)
    get = lambda unit: strings.get((label, unit), [])
    for unit in ("PRO", "PHY", "POT"):
        for closed, pts in get(unit):
            p.path([(q[0], q[1]) for q in pts], closed, fill=ALT_FILL[unit], fill_opacity=0.6, stroke="none")
    for unit, col in LITHO_FILL.items():
        for closed, pts in get(unit):
            p.path([(q[0], q[1]) for q in pts], closed, fill=col, stroke="#444", stroke_width=0.6)
    for unit, col in SHELL_STROKE.items():
        for closed, pts in get(unit):
            p.path([(q[0], q[1]) for q in pts], closed, fill="none", stroke=col, stroke_width=1.8)
    traces = defaultdict(list)
    with open(os.path.join(data, "drillholes", "samples_desurveyed.csv")) as fh:
        for r in csv.DictReader(fh):
            traces[r["HOLE_ID"]].append((float(r["EAST"]), float(r["NORTH"])))
    with open(os.path.join(data, "drillholes", "collar.csv")) as fh:
        for r in csv.DictReader(fh):
            hid = r["HOLE_ID"]
            p.path(traces[hid], False, fill="none", stroke="#333", stroke_width=0.8)
            x, y = p.xy(float(r["EAST"]), float(r["NORTH"]))
            col = "#000" if r["HOLE_TYPE"] == "DD" else "#fff"
            p.parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{col}" stroke="#000" stroke-width="0.8"/>')
    body = p.frame(f"Level plan {label} — drillhole traces, geology, alteration and CuT shells", "Easting (m)", "Northing (m)", 250)
    body += '<text x="868" y="62" font-size="12" font-weight="bold">Holes</text>'
    body += ('<circle cx="875" cy="80" r="3" fill="#000"/><text x="888" y="84" font-size="11">DD</text>'
             '<circle cx="875" cy="97" r="3" fill="#fff" stroke="#000"/><text x="888" y="101" font-size="11">RC</text>')
    body += legend(868, 120, [("QMP", LITHO_FILL["QMP"]), ("IMP", LITHO_FILL["IMP"]), ("BXH", LITHO_FILL["BXH"]),
                              ("DAC", LITHO_FILL["DAC"]), ("Phyllic", ALT_FILL["PHY"]), ("Potassic", ALT_FILL["POT"]),
                              ("Propylitic", ALT_FILL["PRO"])])
    body += legend(868, 250, [("CuT 0.2 %", SHELL_STROKE["CU_0.2"]), ("CuT 0.4 %", SHELL_STROKE["CU_0.4"]),
                              ("CuT 0.8 %", SHELL_STROKE["CU_0.8"])], "line")
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
            f'font-family="Helvetica, Arial, sans-serif"><rect width="100%" height="100%" fill="#fff"/>{body}</svg>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="docs")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    strings = read_strings(os.path.join(a.data, "sections", "interpretation_strings.csv"))
    with open(os.path.join(a.out, "section_+0N.svg"), "w") as fh:
        fh.write(section_svg(a.data, strings, "+0N", N0))
    with open(os.path.join(a.out, "plan_RL1000.svg"), "w") as fh:
        fh.write(plan_svg(a.data, strings, "RL1000"))


if __name__ == "__main__":
    main()
