"""Drillhole planning, downhole surveys, desurveying and interval logging."""
import numpy as np


def direction(az, dip):
    """Unit vector for azimuth (deg from N) and dip (deg, negative = down)."""
    a, d = np.radians(az), np.radians(dip)
    return np.stack([np.sin(a) * np.cos(d), np.cos(a) * np.cos(d), np.sin(d)], axis=-1)


def min_curvature(depths, az, dip, collar):
    """Station coordinates by the minimum-curvature method."""
    t = direction(np.asarray(az), np.asarray(dip))
    xyz = [np.asarray(collar, float)]
    for i in range(1, len(depths)):
        md = depths[i] - depths[i - 1]
        cosb = np.clip(np.dot(t[i - 1], t[i]), -1, 1)
        beta = np.arccos(cosb)
        rf = 1.0 if beta < 1e-9 else 2.0 / beta * np.tan(beta / 2.0)
        xyz.append(xyz[-1] + md / 2.0 * (t[i - 1] + t[i]) * rf)
    return np.array(xyz)


def desurvey(depths, stations_xyz, station_depths):
    out = np.empty((len(depths), 3))
    for k in range(3):
        out[:, k] = np.interp(depths, station_depths, stations_xyz[:, k])
    return out


def plan_holes(model, rng):
    """Drill pattern: E-W sections every 100 m, 100 m spacing in the core,
    200 m step-outs, vertical holes over the centre and -60 angled holes
    drilled in the section plane elsewhere."""
    holes = []
    n_dd = n_rc = 0
    for y in np.arange(-700, 701, 100.0):
        for x in np.arange(-800, 801, 100.0):
            r = np.hypot(x / 1.1, y)
            core = r <= 450
            step = 450 < r <= 820 and (int(x) // 100 + int(y) // 100) % 2 == 0
            if not (core or step):
                continue
            jx, jy = rng.normal(0, 6, 2)
            cx, cy = x + jx, y + jy
            if core and r < 220:
                az, dip = rng.uniform(0, 360), -90.0
            else:
                az = 270.0 if x > 0 else 90.0  # drill toward the stock in the section plane
                az += rng.normal(0, 2)
                dip = -60.0 + rng.normal(0, 1.5)
            is_dd = core or rng.random() < 0.25
            if is_dd:
                n_dd += 1
                hid = f"PDD{n_dd:03d}"
                depth = float(np.round(rng.uniform(520, 820) if core else rng.uniform(380, 560), 1))
            else:
                n_rc += 1
                hid = f"PRC{n_rc:03d}"
                depth = float(np.round(rng.uniform(220, 360) / 2) * 2)
            z = float(model.topo(cx, cy))
            year = 2018 + int(rng.integers(0, 7))
            holes.append(dict(hole_id=hid, x=cx, y=cy, z=z, depth=depth, az=az % 360, dip=dip,
                              type="DD" if is_dd else "RC", year=year, section=f"{int(y):+d}N"))
    return holes


def survey_hole(h, rng):
    """Synthetic downhole survey (gyro) every 30 m with realistic deviation."""
    depths = list(np.arange(0.0, h["depth"], 30.0)) + [h["depth"]]
    depths = np.unique(np.round(depths, 1))
    if h["dip"] <= -89:
        drift_dip = rng.uniform(0.2, 0.6) / 100.0
        drift_az = 0.0
    else:
        drift_dip = rng.uniform(0.8, 2.0) / 100.0  # holes lift (flatten) with depth
        drift_az = rng.normal(0, 1.2) / 100.0
    az = h["az"] + drift_az * depths + rng.normal(0, 0.3, len(depths)).cumsum() * 0.2
    dip = np.minimum(h["dip"] + drift_dip * depths + rng.normal(0, 0.15, len(depths)), -45.0)
    az[0], dip[0] = h["az"], h["dip"]
    return depths, az % 360, dip


def log_intervals(bounds, codes, min_len=0.0):
    """Merge a fine downhole code sequence into from/to intervals.
    bounds has len(codes) + 1 entries (sample boundaries)."""
    out = []
    start = 0
    for i in range(1, len(codes) + 1):
        if i == len(codes) or codes[i] != codes[start]:
            out.append([bounds[start], bounds[i], codes[start]])
            start = i
    if min_len > 0 and len(out) > 1:  # absorb sub-loggable slivers into the previous unit
        merged = [out[0]]
        for fr, to, c in out[1:]:
            if to - fr < min_len:
                merged[-1][1] = to
            elif c == merged[-1][2]:
                merged[-1][1] = to
            else:
                merged.append([fr, to, c])
        out = merged
    return out
