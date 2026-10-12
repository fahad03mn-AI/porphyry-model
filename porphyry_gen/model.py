"""Analytical 3D geological model of a porphyry Cu-Au-Mo deposit.

All coordinates are local metres relative to the stock centre (x = east,
y = north); z is absolute elevation (RL, m).  Every geological unit is
described by an implicit field that is positive inside the unit and roughly
equal to the distance (m) to its boundary, which lets the same model drive
drillhole logging, section strings and wireframe extraction.
"""
import numpy as np

# Projected-coordinate offsets applied when writing files (fictional UTM).
E0 = 352_000.0
N0 = 7_456_000.0

LITHO_CODES = {
    "COL": "Colluvium / post-mineral gravel cover",
    "DAC": "Post-mineral dacite dyke (barren)",
    "BXH": "Hydrothermal breccia pipe (intra-mineral, high grade)",
    "IMP": "Inter-mineral porphyry dyke (diluting, lower grade)",
    "QMP": "Quartz-monzonite porphyry stock (early, main mineralising intrusion)",
    "AND": "Andesite volcanic host rock",
}
LITHO_PRECEDENCE = ["COL", "DAC", "BXH", "IMP", "QMP"]  # youngest first; AND = default

ALT_CODES = {
    "ARG": "Argillic (kaolinite-smectite, near-surface)",
    "PHY": "Phyllic / QSP (quartz-sericite-pyrite)",
    "POT": "Potassic (K-feldspar-biotite-magnetite)",
    "PRO": "Propylitic (chlorite-epidote-calcite)",
    "WK": "Weak / unaltered",
}
ALT_PRECEDENCE = ["ARG", "PHY", "POT", "PRO"]

ZONE_CODES = {
    "COV": "Transported cover (gravels)",
    "LCH": "Leached cap (hematite-goethite-jarosite, Cu depleted)",
    "OXD": "Oxide zone (chrysocolla, brochantite, Cu wad)",
    "SUP": "Supergene enrichment blanket (chalcocite-covellite)",
    "HYP": "Hypogene / primary sulphides (chalcopyrite-bornite)",
}


class SmoothNoise:
    """Deterministic, band-limited random field (sum of random plane waves)."""

    def __init__(self, rng, n=16, wl=(150.0, 600.0), dims=3):
        self.dims = dims
        k = rng.normal(size=(n, 3))
        if dims == 2:
            k[:, 2] = 0.0
        k /= np.linalg.norm(k, axis=1, keepdims=True)
        lam = rng.uniform(*wl, size=n)
        self.k = k * (2 * np.pi / lam)[:, None]
        self.phi = rng.uniform(0, 2 * np.pi, size=n)
        self.scale = np.sqrt(2.0 / n)
        self._cache = {}

    def __call__(self, x, y, z=0.0):
        x, y, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
        # the model re-evaluates the same points many times; cache by content
        key = (x.shape, hash(x.tobytes()), hash(y.tobytes()), 0 if self.dims == 2 else hash(z.tobytes()))
        if key in self._cache:
            return self._cache[key]
        pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
        out = np.empty(len(pts))
        for i in range(0, len(pts), 200_000):
            out[i:i + 200_000] = np.cos(pts[i:i + 200_000] @ self.k.T + self.phi).sum(axis=1)
        out = out.reshape(x.shape) * self.scale
        if len(self._cache) > 6:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = out
        return out


class PorphyryModel:
    STOCK_TOP = 1180.0
    MODEL_ZMIN = 300.0
    MODEL_XY = 1000.0

    def __init__(self, seed=42):
        rng = np.random.default_rng(seed)
        self.n_topo = SmoothNoise(rng, wl=(250, 900), dims=2)
        self.n_cover = SmoothNoise(rng, wl=(150, 500), dims=2)
        self.n_stock = SmoothNoise(rng, wl=(200, 700))
        self.n_imp = SmoothNoise(rng, wl=(150, 500))
        self.n_bx = SmoothNoise(rng, wl=(80, 250))
        self.n_w1 = SmoothNoise(rng, wl=(200, 700), dims=2)
        self.n_w2 = SmoothNoise(rng, wl=(150, 500), dims=2)
        self.n_w3 = SmoothNoise(rng, wl=(200, 700), dims=2)
        self.n_cu = SmoothNoise(rng, wl=(100, 400))
        self.n_mo = SmoothNoise(rng, wl=(100, 400))

    # ------------------------------------------------------------------ surfaces
    def topo(self, x, y):
        hill = 70.0 * np.exp(-((x - 150) ** 2 + (y + 100) ** 2) / (2 * 550.0 ** 2))
        return 1240.0 + hill + 0.02 * x - 0.015 * y + 9.0 * self.n_topo(x, y)

    def cover_thickness(self, x, y):
        # thicker gravels on the lower (north-west) flank
        base = 6.0 + 8.0 * np.clip((-(x) + y) / 1400.0 + 0.3, 0, 1)
        return np.clip(base + 4.0 * self.n_cover(x, y), 2.0, None)

    def _phyllic_surface_factor(self, x, y):
        rho = self.stock_rho(x, y, self.topo(x, y) - 60.0)
        return np.exp(-(rho / 2.0) ** 2)

    def weathering_depths(self, x, y):
        """Depth below topography of BOLC, BOX and TOH (m)."""
        pf = self._phyllic_surface_factor(x, y)
        d_bolc = 35.0 + 55.0 * pf + 12.0 * self.n_w1(x, y)
        d_bolc = np.maximum(d_bolc, self.cover_thickness(x, y) + 8.0)
        d_box = d_bolc + np.maximum(22.0 + 8.0 * self.n_w2(x, y), 6.0)  # oxide >= 6 m
        d_toh = d_box + np.maximum(30.0 + 75.0 * pf + 12.0 * self.n_w3(x, y), 10.0)  # blanket >= 10 m
        return d_bolc, d_box, d_toh, pf

    # ------------------------------------------------------------------ stock geometry
    def stock_axis(self, z):
        h = self.STOCK_TOP - z
        return 0.08 * h, 0.05 * h  # axis plunges gently to the NE with depth

    def stock_axes(self, z):
        h = np.clip(self.STOCK_TOP - z, 0, None)
        a = 190.0 + 0.25 * h
        return a, 0.72 * a

    def stock_rho(self, x, y, z):
        """Normalised elliptical radius from the stock axis (1 = stock wall)."""
        cx, cy = self.stock_axis(z)
        a, b = self.stock_axes(z)
        th = np.radians(30.0)  # long axis trends 060
        dx, dy = x - cx, y - cy
        u = dx * np.cos(th) + dy * np.sin(th)
        v = -dx * np.sin(th) + dy * np.cos(th)
        rho = np.sqrt((u / a) ** 2 + (v / b) ** 2)
        return rho * (1.0 + 0.07 * self.n_stock(x, y, z))

    def stock_rscale(self, z):
        a, b = self.stock_axes(z)
        return np.sqrt(a * b)

    # ------------------------------------------------------------------ unit fields (+ inside)
    def f_qmp(self, x, y, z):
        rho = self.stock_rho(x, y, z)
        cap = self.STOCK_TOP - 120.0 * rho ** 2 - z
        return np.minimum((1.0 - rho) * self.stock_rscale(z), cap)

    # strike (deg), centre point and along-strike half-length used for aligned meshing
    DYKES = {"IMP": (50.0, (40.0, 20.0), 700.0, (300.0, 1160.0)),
             "DAC": (135.0, (-100.0, 150.0), 860.0, (300.0, 1420.0))}

    def f_imp(self, x, y, z):
        # dyke striking 050, dipping 80 SE, rooted in the stock
        strike, dip = np.radians(50.0), np.radians(80.0)
        n = np.array([np.cos(strike) * np.sin(dip), -np.sin(strike) * np.sin(dip), np.cos(dip)])
        d = (x - 40.0) * n[0] + (y - 20.0) * n[1] + (z - 1000.0) * n[2]
        hw = 18.0 + 5.0 * self.n_imp(x, y, z)
        rho = self.stock_rho(x, y, z)
        return np.minimum.reduce([hw - np.abs(d), 1150.0 - z, (2.0 - rho) * 120.0])

    def f_dac(self, x, y, z):
        # post-mineral dyke striking 135, dipping 85 NE, cuts everything below cover
        strike, dip = np.radians(135.0), np.radians(85.0)
        n = np.array([np.cos(strike) * np.sin(dip), -np.sin(strike) * np.sin(dip), np.cos(dip)])
        sx, sy = np.sin(strike), np.cos(strike)
        d = (x + 100.0) * n[0] + (y - 150.0) * n[1] + (z - 1000.0) * n[2]
        along = (x + 100.0) * sx + (y - 150.0) * sy
        top = self.topo(x, y) - self.cover_thickness(x, y) - z
        return np.minimum.reduce([10.0 - np.abs(d), 850.0 - np.abs(along), top])

    def f_bxh(self, x, y, z):
        cx, cy = 120.0 + 0.10 * (1230.0 - z), -80.0 - 0.05 * (1230.0 - z)
        r = np.hypot(x - cx, y - cy)
        rad = 45.0 + 10.0 * self.n_bx(x, y, z) + 0.06 * np.clip(z - 900.0, 0, None)
        return np.minimum.reduce([rad - r, z - 900.0, 1230.0 - z])

    def f_col(self, x, y, z):
        t = self.topo(x, y)
        return np.minimum(self.cover_thickness(x, y) - (t - z), t - z)

    def litho_fields(self, x, y, z):
        return {
            "COL": self.f_col(x, y, z),
            "DAC": self.f_dac(x, y, z),
            "BXH": self.f_bxh(x, y, z),
            "IMP": self.f_imp(x, y, z),
            "QMP": self.f_qmp(x, y, z),
        }

    def f_pot(self, x, y, z):
        rho = self.stock_rho(x, y, z)
        top = 1150.0 - 40.0 * rho + 30.0 * self.n_imp(x, y, 0.0 * z)
        return np.minimum((1.5 - rho) * self.stock_rscale(z), top - z)

    def f_phy(self, x, y, z):
        rho = self.stock_rho(x, y, z)
        zb = 1080.0 + 35.0 * self.n_bx(x, y, 0.0 * z) - 190.0 * np.clip(rho - 0.8, 0, None)
        return np.minimum((2.4 - rho) * self.stock_rscale(z), z - zb)

    def f_arg(self, x, y, z):
        rho = self.stock_rho(x, y, z)
        return np.minimum(70.0 - (self.topo(x, y) - z), (1.9 - rho) * self.stock_rscale(z))

    def f_pro(self, x, y, z):
        return (4.2 - self.stock_rho(x, y, z)) * self.stock_rscale(z)

    def alt_fields(self, x, y, z):
        return {"ARG": self.f_arg(x, y, z), "PHY": self.f_phy(x, y, z),
                "POT": self.f_pot(x, y, z), "PRO": self.f_pro(x, y, z)}

    # ------------------------------------------------------------------ classification
    @staticmethod
    def _classify(fields, order, default):
        shape = next(iter(fields.values())).shape
        out = np.full(shape, default, dtype=object)
        for code in reversed(order):  # oldest first, younger overwrite
            out[fields[code] > 0] = code
        return out

    def lithology(self, x, y, z):
        return self._classify(self.litho_fields(x, y, z), LITHO_PRECEDENCE, "AND")

    def alteration(self, x, y, z):
        alt = self._classify(self.alt_fields(x, y, z), ALT_PRECEDENCE, "WK")
        lit = self.lithology(x, y, z)
        alt[(lit == "DAC") | (lit == "COL")] = "WK"  # post-mineral units are fresh
        return alt

    def zone(self, x, y, z):
        depth = self.topo(x, y) - z
        d_bolc, d_box, d_toh, _ = self.weathering_depths(x, y)
        out = np.full(np.shape(depth), "HYP", dtype=object)
        out[depth < d_toh] = "SUP"
        out[depth < d_box] = "OXD"
        out[depth < d_bolc] = "LCH"
        out[self.f_col(x, y, z) > 0] = "COV"
        return out

    # ------------------------------------------------------------------ grades (expected values)
    def grades(self, x, y, z):
        """Expected (noise-free at sample scale) grades.  Returns dict of arrays."""
        x, y, z = np.broadcast_arrays(*(np.asarray(v, float) for v in (x, y, z)))
        rho = self.stock_rho(x, y, z)
        rs = self.stock_rscale(z)
        dm = (rho - 1.0) * rs  # m outside the stock wall (negative inside)
        vf = np.exp(-((z - 1080.0) / 380.0) ** 2)
        width = np.where(dm < 20.0, 150.0, 115.0)
        shell = 0.85 * np.exp(-((dm - 20.0) / width) ** 2)
        core = 1.0 - 0.6 * np.clip((-dm - 120.0) / 150.0, 0, 1) * np.clip((1050.0 - z) / 200.0, 0, 1)
        cu_h = (0.03 + shell * vf * core) * np.exp(0.30 * self.n_cu(x, y, z))

        au_h = cu_h * 0.42 * np.where(dm < 0, 1.3, 0.9)
        mo = (15.0 + 260.0 * np.exp(-((dm - 130.0) / 130.0) ** 2) * np.exp(-((z - 850.0) / 320.0) ** 2)) \
            * np.exp(0.35 * self.n_mo(x, y, z))
        ag = 0.4 + 2.6 * cu_h

        lit = self.lithology(x, y, z)
        fac = {"IMP": (0.45, 0.4, 0.6), "BXH": (1.6, 1.8, 1.2), "DAC": (0.0, 0.0, 0.0), "COL": (0.0, 0.0, 0.0)}
        for code, (fc, fa, fm) in fac.items():
            m = lit == code
            cu_h[m] *= fc
            au_h[m] *= fa
            mo[m] *= fm
            ag[m] *= fc if fc > 0 else 1.0
        post = (lit == "DAC") | (lit == "COL")
        cu_h[post] = np.where(lit[post] == "COL", 0.02, 0.006)
        au_h[post] = 0.008
        mo[post] = 2.0
        ag[post] = 0.3

        alt = self.alteration(x, y, z)
        s_by_alt = {"PHY": 4.5, "ARG": 2.2, "POT": 1.6, "PRO": 1.1, "WK": 0.3}
        s = np.zeros(z.shape)
        for code, v in s_by_alt.items():
            s[alt == code] = v

        # ---- supergene profile
        depth = self.topo(x, y) - z
        d_bolc, d_box, d_toh, pf = self.weathering_depths(x, y)
        zone = self.zone(x, y, z)
        cu = cu_h.copy()
        au = au_h.copy()
        cu_as = 0.03 * cu
        cu_cn = 0.06 * cu
        m = zone == "LCH"
        cu[m] *= 0.08
        au[m] *= 0.95
        mo[m] *= 0.6
        s[m] = 0.15
        cu_as[m], cu_cn[m] = 0.30 * cu[m], 0.15 * cu[m]
        m = zone == "OXD"
        cu[m] *= 0.95
        s[m] = 0.4
        cu_as[m], cu_cn[m] = 0.72 * cu[m], 0.10 * cu[m]
        m = zone == "SUP"
        frac = np.clip((depth - d_box) / np.maximum(d_toh - d_box, 1.0), 0, 1)
        enr = 1.0 + 1.7 * (1.0 - frac) ** 1.5 * (0.6 + 0.4 * pf)
        cu[m] *= enr[m]
        cu_as[m], cu_cn[m] = 0.05 * cu[m], (0.70 - 0.35 * frac[m]) * cu[m]
        m = zone == "COV"
        cu_as[m], cu_cn[m] = 0.2 * cu[m], 0.1 * cu[m]

        return {"CuT_pct": cu, "CuAS_pct": cu_as, "CuCN_pct": cu_cn, "Au_gpt": au,
                "Mo_ppm": mo, "Ag_gpt": ag, "S_pct": s, "_lit": lit, "_alt": alt, "_zone": zone}

    def density(self, lit, zone):
        sg = np.full(np.shape(lit), 2.68)
        for code, v in {"QMP": 2.62, "IMP": 2.60, "BXH": 2.55, "DAC": 2.58, "AND": 2.70}.items():
            sg[lit == code] = v
        for code, v in {"LCH": 2.25, "OXD": 2.40, "SUP": 2.52}.items():
            sg[zone == code] = v
        sg[lit == "COL"] = 1.95
        return sg
