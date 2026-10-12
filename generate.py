#!/usr/bin/env python3
"""Generate a synthetic porphyry Cu-Au-Mo dataset for wireframe-modelling practice.

    python3 generate.py                 # default seed 42, writes ./data
    python3 generate.py --seed 7 --out data_seed7 --dxf-solids

Only numpy is required.
"""
import argparse
import os

import numpy as np

from porphyry_gen import drillholes as dh
from porphyry_gen.meshing import (grid_surface, is_closed, marching_squares, marching_tetrahedra,
                                  mesh_area, mesh_volume)
from porphyry_gen.model import (ALT_CODES, E0, LITHO_CODES, N0, ZONE_CODES, PorphyryModel)
from porphyry_gen.writers import (write_csv, write_dxf_mesh, write_dxf_polylines, write_obj,
                                  write_surpac_str)

STR_NUMBERS = {"TOPO": 1, "BOC": 2, "BOLC": 3, "BOX": 4, "TOH": 5,
               "QMP": 10, "IMP": 11, "BXH": 12, "DAC": 13,
               "POT": 20, "PHY": 21, "ARG": 22, "PRO": 23,
               "CU_0.2": 30, "CU_0.4": 31, "CU_0.8": 32}

ASSAY_FIELDS = [("CuT_pct", 3, 0.001, 0.22), ("CuAS_pct", 3, 0.001, 0.25), ("CuCN_pct", 3, 0.001, 0.25),
                ("Au_gpt", 3, 0.005, 0.45), ("Mo_ppm", 0, 1, 0.40), ("Ag_gpt", 1, 0.2, 0.35),
                ("S_pct", 2, 0.01, 0.20)]


def offset(xyz):
    out = np.array(xyz, float, copy=True)
    out[..., 0] += E0
    out[..., 1] += N0
    return out


# ---------------------------------------------------------------------------- drillholes
def build_drillholes(m, rng, out):
    d = os.path.join(out, "drillholes")
    os.makedirs(d, exist_ok=True)
    holes = dh.plan_holes(m, rng)
    collar, survey, assay, litho, alter, zone, dens, contacts, desurv = ([] for _ in range(9))
    sample_no = 100000
    for h in holes:
        collar.append([h["hole_id"], f"{h['x'] + E0:.2f}", f"{h['y'] + N0:.2f}", f"{h['z']:.2f}",
                       f"{h['depth']:.1f}", f"{h['az']:.1f}", f"{h['dip']:.1f}", h["type"], h["year"], h["section"]])
        sd, saz, sdip = dh.survey_hole(h, rng)
        survey += [[h["hole_id"], f"{a:.1f}", f"{b:.1f}", f"{c:.1f}"] for a, b, c in zip(sd, saz, sdip)]
        st_xyz = dh.min_curvature(sd, saz, sdip, (h["x"], h["y"], h["z"]))

        # --- geological logging at 0.5 m resolution
        b = np.round(np.append(np.arange(0.0, h["depth"], 0.5), h["depth"]), 2)
        mid = (b[:-1] + b[1:]) / 2.0
        p = dh.desurvey(mid, st_xyz, sd)
        lit = m.lithology(*p.T)
        alt = m.alteration(*p.T)
        zon = m.zone(*p.T)
        afields = m.alt_fields(*p.T)
        for fr, to, c in dh.log_intervals(b, lit, min_len=1.0):
            litho.append([h["hole_id"], f"{fr:.1f}", f"{to:.1f}", c])
        for fr, to, c in dh.log_intervals(b, alt, min_len=1.0):
            sel = (mid >= fr) & (mid < to)
            if c == "WK":
                inten = 0
            else:  # intensity from distance inside the alteration envelope
                depth_in = np.median(afields[c][sel])
                inten = 3 if depth_in > 80 else 2 if depth_in > 25 else 1
            alter.append([h["hole_id"], f"{fr:.1f}", f"{to:.1f}", c, inten])
        zint = dh.log_intervals(b, zon)
        for fr, to, c in zint:
            zone.append([h["hole_id"], f"{fr:.1f}", f"{to:.1f}", c])

        # --- contact points (for implicit modelling / snapping section strings)
        surf_name = {("COV", "LCH"): "BOC", ("LCH", "OXD"): "BOLC", ("OXD", "SUP"): "BOX", ("SUP", "HYP"): "TOH"}
        lints = dh.log_intervals(b, lit, min_len=1.0)
        for kind, ints in (("LITHOLOGY", lints), ("WEATHERING", zint)):
            for (f0, t0, c0), (f1, t1, c1) in zip(ints[:-1], ints[1:]):
                x, y, z = offset(dh.desurvey([t0], st_xyz, sd)[0])
                name = surf_name.get((c0, c1), f"{c0}-{c1}") if kind == "WEATHERING" else f"{c0}-{c1}"
                contacts.append([h["hole_id"], f"{t0:.1f}", f"{x:.2f}", f"{y:.2f}", f"{z:.2f}", kind, name, c0, c1])

        # --- 2 m assay samples, composited from 4 sub-points for sample support
        sb = np.append(np.arange(0.0, h["depth"], 2.0), h["depth"])
        sb = sb[np.append(np.diff(sb) > 0.05, True)]
        fr_, to_ = sb[:-1], sb[1:]
        ns = len(fr_)
        sub = fr_[:, None] + (to_ - fr_)[:, None] * (np.arange(4) + 0.5) / 4.0
        g = m.grades(*dh.desurvey(sub.ravel(), st_xyz, sd).T)
        gm = {k: g[k].reshape(ns, 4).mean(axis=1) for k, *_ in ASSAY_FIELDS}
        pm = dh.desurvey((fr_ + to_) / 2.0, st_xyz, sd)
        lit_m, alt_m, zon_m = m.lithology(*pm.T), m.alteration(*pm.T), m.zone(*pm.T)
        vals = np.stack([gm[k] * np.exp(rng.normal(0, sig, ns)) for k, _, _, sig in ASSAY_FIELDS], axis=1)
        sol = vals[:, 1] + vals[:, 2]
        scale = np.where(sol > 0.98 * vals[:, 0], 0.98 * vals[:, 0] / sol, 1.0)
        vals[:, 1] *= scale
        vals[:, 2] *= scale
        vals = np.maximum(vals, np.array([dl / 2.0 for _, _, dl, _ in ASSAY_FIELDS]))
        unsampled = ((h["type"] == "RC") & (lit_m == "COL")) | (rng.random(ns) < 0.004)
        rec = np.clip(rng.normal(np.where(zon_m == "HYP", 97, 88), 4), 55, 100)
        for i in range(ns):
            sample_no += 1
            sid = f"S{sample_no}"
            fmt = ["" if unsampled[i] else (f"{v:.{dp}f}" if dp else f"{round(v):d}")
                   for v, (_, dp, _, _) in zip(vals[i], ASSAY_FIELDS)]
            r = "" if h["type"] == "RC" else f"{rec[i]:.0f}"
            assay.append([h["hole_id"], sid, f"{fr_[i]:.1f}", f"{to_[i]:.1f}", r] + fmt)
            x, y, z = offset(pm[i])
            desurv.append([h["hole_id"], sid, f"{fr_[i]:.1f}", f"{to_[i]:.1f}", f"{x:.2f}", f"{y:.2f}", f"{z:.2f}",
                           lit_m[i], alt_m[i], zon_m[i]] + fmt)

        # --- density (DD only, ~every 20 m)
        if h["type"] == "DD":
            dd = np.arange(rng.uniform(5, 20), h["depth"], 20.0)
            pp = dh.desurvey(dd, st_xyz, sd)
            lc, zc = m.lithology(*pp.T), m.zone(*pp.T)
            sg = m.density(lc, zc) + rng.normal(0, 0.035, len(dd))
            for a_, b_, c_ in zip(dd, sg, zc):
                dens.append([h["hole_id"], f"{a_:.2f}", f"{b_:.3f}",
                             "Archimedes (wax-coated)" if c_ in ("COV", "LCH", "OXD") else "Archimedes"])

    write_csv(f"{d}/collar.csv", ["HOLE_ID", "EAST", "NORTH", "RL", "DEPTH", "AZIMUTH", "DIP", "HOLE_TYPE", "YEAR", "SECTION"], collar)
    write_csv(f"{d}/survey.csv", ["HOLE_ID", "DEPTH", "AZIMUTH", "DIP"], survey)
    write_csv(f"{d}/assay.csv", ["HOLE_ID", "SAMPLE_ID", "FROM", "TO", "RECOVERY_PCT"] + [a[0] for a in ASSAY_FIELDS], assay)
    write_csv(f"{d}/lithology.csv", ["HOLE_ID", "FROM", "TO", "LITHO"], litho)
    write_csv(f"{d}/alteration.csv", ["HOLE_ID", "FROM", "TO", "ALTERATION", "INTENSITY"], alter)
    write_csv(f"{d}/weathering.csv", ["HOLE_ID", "FROM", "TO", "MIN_ZONE"], zone)
    write_csv(f"{d}/density.csv", ["HOLE_ID", "DEPTH", "SG", "METHOD"], dens)
    write_csv(f"{d}/contacts.csv", ["HOLE_ID", "DEPTH", "EAST", "NORTH", "RL", "CONTACT_TYPE", "SURFACE", "UNIT_ABOVE", "UNIT_BELOW"], contacts)
    write_csv(f"{d}/samples_desurveyed.csv", ["HOLE_ID", "SAMPLE_ID", "FROM", "TO", "EAST", "NORTH", "RL", "LITHO", "ALTERATION", "MIN_ZONE"]
              + [a[0] for a in ASSAY_FIELDS], desurv)
    return holes, len(assay)


# ---------------------------------------------------------------------------- wireframes
def grid3(m, cell, bbox):
    (x0, x1), (y0, y1), (z0, z1) = bbox
    xs, ys, zs = np.arange(x0, x1 + 1e-6, cell), np.arange(y0, y1 + 1e-6, cell), np.arange(z0, z1 + 1e-6, cell)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return xs, ys, zs, X, Y, Z


def build_wireframes(m, out, dxf_solids=False):
    sd = os.path.join(out, "wireframes", "surfaces")
    so = os.path.join(out, "wireframes", "solids")
    os.makedirs(sd, exist_ok=True)
    os.makedirs(so, exist_ok=True)
    report = []

    # --- DTM surfaces
    xs = ys = np.arange(-1000.0, 1000.0 + 1e-6, 25.0)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    t = m.topo(X, Y)
    d_bolc, d_box, d_toh, _ = m.weathering_depths(X, Y)
    surfaces = {"topography": t, "base_of_cover_BOC": t - m.cover_thickness(X, Y),
                "base_of_leached_cap_BOLC": t - d_bolc, "base_of_oxide_BOX": t - d_box,
                "top_of_hypogene_TOH": t - d_toh}
    for name, zg in surfaces.items():
        v, tr = grid_surface(xs, ys, zg)
        v = offset(v)
        write_obj(f"{sd}/{name}.obj", v, tr, name)
        write_dxf_mesh(f"{sd}/{name}.dxf", v, tr, name.upper())
        report.append([name, "surface", len(tr), f"{mesh_area(v, tr):.0f}", "", "", "", ""])

    # --- solids on a 20 m grid: (name, field, cell, origin, rotation)
    full = ((-1000, 1000), (-1000, 1000), (300, 1420))
    _, _, _, X, Y, Z = grid3(m, 20.0, full)
    org = (full[0][0], full[1][0], full[2][0])
    topo_clip = m.topo(X, Y) - Z
    gr = m.grades(X, Y, Z)
    solids = [
        ("litho_QMP_stock", np.minimum(m.f_qmp(X, Y, Z), topo_clip), 20.0, org, None),
        ("alt_POT_potassic", np.minimum(m.f_pot(X, Y, Z), topo_clip), 20.0, org, None),
        ("alt_PHY_phyllic", np.minimum(m.f_phy(X, Y, Z), topo_clip), 20.0, org, None),
        ("alt_ARG_argillic", np.minimum(m.f_arg(X, Y, Z), topo_clip), 20.0, org, None),
        ("alt_PRO_propylitic_outer", np.minimum(m.f_pro(X, Y, Z), topo_clip)[::2, ::2, ::2], 40.0, org, None),
    ]
    for cut in (0.2, 0.4, 0.8):
        solids.append((f"gradeshell_CuT_{cut:.1f}pct", np.minimum(gr["CuT_pct"] - cut, topo_clip / 100.0), 20.0, org, None))
    solids.append(("gradeshell_Au_0.3gpt", np.minimum(gr["Au_gpt"] - 0.3, topo_clip / 100.0), 20.0, org, None))
    solids.append(("gradeshell_Mo_100ppm", np.minimum(gr["Mo_ppm"] - 100.0, topo_clip)[::2, ::2, ::2], 40.0, org, None))

    # reference tonnes / grade inside each grade shell (20 m point support)
    tonnes = m.density(gr["_lit"], gr["_zone"]) * 20.0 ** 3 * (topo_clip > 0)
    ref = {}
    for name, el, cut in [(s[0], *k) for s, k in zip(solids[5:], [("CuT_pct", 0.2), ("CuT_pct", 0.4), ("CuT_pct", 0.8),
                                                                  ("Au_gpt", 0.3), ("Mo_ppm", 100.0)])]:
        t = tonnes * (gr[el] >= cut)
        tn = t.sum()
        ref[name] = (f"{tn / 1e6:.1f}", f"{(t * gr['CuT_pct']).sum() / tn:.3f} / {(t * gr['Au_gpt']).sum() / tn:.3f} / "
                                        f"{(t * gr['Mo_ppm']).sum() / tn:.0f}")
    del gr, X, Y, Z, topo_clip

    # breccia pipe: small 8 m grid
    bb = ((-100, 300), (-300, 100), (880, 1250))
    _, _, _, Xf, Yf, Zf = grid3(m, 8.0, bb)
    solids.append(("litho_BXH_breccia_pipe", np.minimum(m.f_bxh(Xf, Yf, Zf), m.topo(Xf, Yf) - Zf), (8.0, 8.0, 8.0),
                   (bb[0][0], bb[1][0], bb[2][0]), None))
    # thin dykes: grid aligned with strike (25 m along, 4 m across, 20 m vertical)
    for code, fn in (("IMP", m.f_imp), ("DAC", m.f_dac)):
        strike, (cx, cy), half, (z0, z1) = m.DYKES[code]
        s_ = np.radians(strike)
        u = np.array([np.sin(s_), np.cos(s_), 0.0])  # along strike
        w = np.array([np.cos(s_), -np.sin(s_), 0.0])  # across strike
        rot = np.stack([u, w, [0.0, 0.0, 1.0]], axis=1)
        us, ws, zs = np.arange(-half - 50, half + 51, 25.0), np.arange(-60.0, 60.1, 4.0), np.arange(z0, z1 + 1e-6, 20.0)
        U, W, Zl = np.meshgrid(us, ws, zs, indexing="ij")
        Xd, Yd = cx + U * u[0] + W * w[0], cy + U * u[1] + W * w[1]
        f = np.minimum(fn(Xd, Yd, Zl), m.topo(Xd, Yd) - Zl)
        name = "litho_IMP_dyke" if code == "IMP" else "litho_DAC_dyke"
        solids.append((name, f, (25.0, 4.0, 20.0), (us[0], ws[0], zs[0]), (rot, (cx, cy, 0.0))))

    for name, f, cell, origin, rot in solids:
        cell = cell if isinstance(cell, tuple) else (cell, cell, cell)
        v, tr = marching_tetrahedra(f, origin, cell)
        if rot is not None:
            R, shift = rot
            v = v @ R.T + np.asarray(shift)
            if np.linalg.det(R) < 0:
                tr = tr[:, ::-1]
        v = offset(v)
        vol = mesh_volume(v, tr)
        write_obj(f"{so}/{name}.obj", v, tr, name)
        if dxf_solids:
            write_dxf_mesh(f"{so}/{name}.dxf", v, tr, name.upper())
        report.append([name, "solid", len(tr), f"{mesh_area(v, tr):.0f}", f"{vol:.0f}",
                       "yes" if is_closed(tr) else "no", *ref.get(name, ("", ""))])
        print(f"  {name:32s} {len(tr):7d} tris  vol {vol / 1e6:9.2f} Mm3")
    return report


# ---------------------------------------------------------------------------- sections & plans
def contour_strings(f2d, coords_fn, min_pts=4):
    out = []
    for pts, closed in marching_squares(f2d):
        if len(pts) >= min_pts:
            out.append((coords_fn(pts), closed))
    return out


def build_sections(m, out):
    d = os.path.join(out, "sections")
    os.makedirs(d, exist_ok=True)
    step = 10.0
    xs = np.arange(-1000.0, 1000.0 + 1e-6, step)
    zs = np.arange(300.0, 1420.0 + 1e-6, step)
    all_strings, csv_rows = [], []
    sid = 0

    def add(kind, unit, label, pts3, closed, section):
        nonlocal sid
        sid += 1
        p = offset(pts3)
        all_strings.append((kind, unit, label, p, closed))
        for n, (x, y, z) in enumerate(p):
            csv_rows.append([sid, section, kind, unit, int(closed), n + 1, f"{x:.2f}", f"{y:.2f}", f"{z:.2f}"])

    def planes(kind, fixed_vals, grid_fn):
        for val in fixed_vals:
            label = f"{int(val):+d}N" if kind == "SECTION" else f"RL{int(val)}"
            X, Y, Z, to3 = grid_fn(val)
            tclip = m.topo(X, Y) - Z
            if kind == "SECTION":  # open strings for topography and weathering surfaces
                y = np.full_like(xs, val)
                t = m.topo(xs, y)
                db, dx_, dt, _ = m.weathering_depths(xs, y)
                for unit, zz in (("TOPO", t), ("BOC", t - m.cover_thickness(xs, y)), ("BOLC", t - db),
                                 ("BOX", t - dx_), ("TOH", t - dt)):
                    add(kind, unit, label, np.stack([xs, y, zz], 1), False, label)
            fields = {"QMP": m.f_qmp(X, Y, Z), "IMP": m.f_imp(X, Y, Z), "BXH": m.f_bxh(X, Y, Z),
                      "DAC": m.f_dac(X, Y, Z), "POT": m.f_pot(X, Y, Z), "PHY": m.f_phy(X, Y, Z),
                      "ARG": m.f_arg(X, Y, Z), "PRO": m.f_pro(X, Y, Z)}
            cu = m.grades(X, Y, Z)["CuT_pct"]
            for c in (0.2, 0.4, 0.8):
                fields[f"CU_{c:.1f}"] = (cu - c) * 100.0
            for unit, f in fields.items():
                for pts3, closed in contour_strings(np.minimum(f, tclip), to3):
                    add(kind, unit, label, pts3, closed, label)

    def section_grid(yv):
        X, Z = np.meshgrid(xs, zs, indexing="ij")
        Y = np.full_like(X, yv)
        to3 = lambda p: np.stack([xs[0] + p[:, 0] * step, np.full(len(p), yv), zs[0] + p[:, 1] * step], 1)
        return X, Y, Z, to3

    pxy = np.arange(-1000.0, 1000.0 + 1e-6, step)

    def plan_grid(zv):
        X, Y = np.meshgrid(pxy, pxy, indexing="ij")
        Z = np.full_like(X, zv)
        to3 = lambda p: np.stack([pxy[0] + p[:, 0] * step, pxy[0] + p[:, 1] * step, np.full(len(p), zv)], 1)
        return X, Y, Z, to3

    planes("SECTION", np.arange(-700.0, 701.0, 100.0), section_grid)
    planes("PLAN", np.arange(500.0, 1201.0, 100.0), plan_grid)

    write_csv(f"{d}/interpretation_strings.csv",
              ["STRING_ID", "SECTION_OR_PLAN", "VIEW", "UNIT", "CLOSED", "SEQ", "EAST", "NORTH", "RL"], csv_rows)
    for kind in ("SECTION", "PLAN"):
        sel = [s for s in all_strings if s[0] == kind]
        tag = "sections_EW" if kind == "SECTION" else "level_plans"
        write_dxf_polylines(f"{d}/{tag}.dxf", [(f"{u}_{lab}".replace(".", "p").replace("+", "P").replace("-", "M"), p, c)
                                               for _, u, lab, p, c in sel])
        groups = {"geology": ("TOPO", "BOC", "BOLC", "BOX", "TOH", "QMP", "IMP", "BXH", "DAC"),
                  "alteration": ("POT", "PHY", "ARG", "PRO"), "cu_gradeshells": ("CU_0.2", "CU_0.4", "CU_0.8")}
        for gname, units in groups.items():
            write_surpac_str(f"{d}/{tag}_{gname}.str", f"{tag}_{gname}",
                             [(STR_NUMBERS[u], p, c, f"{u} {lab}") for _, u, lab, p, c in sel if u in units])
    return len(all_strings)


def write_codes(out):
    rows = [["LITHO", k, v] for k, v in LITHO_CODES.items()]
    rows += [["ALTERATION", k, v] for k, v in ALT_CODES.items()]
    rows += [["MIN_ZONE", k, v] for k, v in ZONE_CODES.items()]
    rows += [["SURPAC_STRING", str(v), k] for k, v in STR_NUMBERS.items()]
    write_csv(os.path.join(out, "codes.csv"), ["TABLE", "CODE", "DESCRIPTION"], rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="data")
    ap.add_argument("--dxf-solids", action="store_true",
                    help="also write solids as DXF 3DFACE (large files; OBJ is always written)")
    a = ap.parse_args()
    m = PorphyryModel(a.seed)
    rng = np.random.default_rng(a.seed + 1)
    os.makedirs(a.out, exist_ok=True)
    print("drillholes ...")
    holes, n_samples = build_drillholes(m, rng, a.out)
    print(f"  {len(holes)} holes, {n_samples} samples")
    print("sections & plans ...")
    print(f"  {build_sections(m, a.out)} strings")
    print("wireframes ...")
    report = build_wireframes(m, a.out, dxf_solids=a.dxf_solids)
    write_csv(os.path.join(a.out, "wireframes", "wireframe_report.csv"),
              ["NAME", "TYPE", "TRIANGLES", "AREA_M2", "VOLUME_M3", "CLOSED", "TRUE_TONNES_Mt", "TRUE_CuT%/Au_gpt/Mo_ppm"], report)
    write_codes(a.out)
    print("done")


if __name__ == "__main__":
    main()
