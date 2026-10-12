"""Iso-surface (marching tetrahedra) and iso-line (marching squares) extraction.

Both work on fields that are positive inside a unit.  Vertices are shared via
grid-edge keys so the resulting meshes are watertight and the polylines are
continuous.
"""
import numpy as np

# Cube corner offsets (i, j, k) and a 6-tetrahedra split along the 0-6 diagonal.
_CORNERS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
_TETS = [(0, 5, 1, 6), (0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6)]


def marching_tetrahedra(f, origin, spacing, iso=0.0, snap=0.3, max_iter=6):
    """Return (vertices[n,3], triangles[m,3]) of the closed surface f = iso.

    f is a (nx, ny, nz) array sampled at origin + index * spacing.  The field is
    padded with an 'outside' value so the surface is always closed.  Crossings
    closer than `snap` (fraction of an edge) to a grid node are snapped onto the
    node, which removes most of the sliver triangles marching tetrahedra makes.
    """
    nosnap = np.zeros(0, dtype=np.int64)
    for _ in range(max_iter):
        verts, tris, keys, n_nodes = _mt(f, origin, spacing, iso, snap, nosnap)
        bad = _nonmanifold_vertices(tris)
        bad_nodes = keys[bad][keys[bad] < n_nodes]  # snapped vertices on pinched edges
        if len(bad_nodes) == 0:
            break
        nosnap = np.union1d(nosnap, bad_nodes)
    return verts, tris


def _nonmanifold_vertices(tris):
    e = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return np.unique(u[c != 2])


def _mt(f, origin, spacing, iso, snap, nosnap):
    f = np.pad(np.asarray(f, float) - iso, 1, constant_values=-1.0)
    f[f == 0.0] = 1e-9
    origin = np.asarray(origin, float) - np.asarray(spacing, float)
    spacing = np.asarray(spacing, float)
    nx, ny, nz = f.shape
    ii, jj, kk = np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), np.arange(nz - 1), indexing="ij")
    base = np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1)
    inside = f > 0
    # keep only cubes that straddle the surface
    cube_in = np.stack([inside[tuple((base + c).T)] for c in _CORNERS], axis=1)
    mixed = cube_in.any(axis=1) & ~cube_in.all(axis=1)
    base = base[mixed]

    def lin(idx):
        return (idx[:, 0] * ny + idx[:, 1]) * nz + idx[:, 2]

    edge_key, edge_pts, tris_keys, tris_inside = [], [], [], []
    for tet in _TETS:
        idx = [base + _CORNERS[c] for c in tet]
        val = np.stack([f[tuple(p.T)] for p in idx], axis=1)
        ins = val > 0
        code = (ins * (1 << np.arange(4))).sum(axis=1)
        ids = np.stack([lin(p) for p in idx], axis=1)
        pos = np.stack(idx, axis=1).astype(float)

        def edge(a, b, sel):
            fa, fb = val[sel, a], val[sel, b]
            t = fa / (fa - fb)
            ka, kb = ids[sel, a], ids[sel, b]
            key = np.minimum(ka, kb) * (nx * ny * nz) + np.maximum(ka, kb) + nx * ny * nz
            at_a = (t < snap) & ~np.isin(ka, nosnap)
            at_b = (t > 1.0 - snap) & ~np.isin(kb, nosnap)
            t = np.where(at_a, 0.0, np.where(at_b, 1.0, t))[:, None]
            key = np.where(at_a, ka, np.where(at_b, kb, key))  # snapped: key = node id
            p = pos[sel, a] + t * (pos[sel, b] - pos[sel, a])
            return key, p

        for c in range(1, 15):
            sel = code == c
            if not sel.any():
                continue
            in_c = [v for v in range(4) if c >> v & 1]
            out_c = [v for v in range(4) if not c >> v & 1]
            # direction pointing from the inside corners to the outside corners
            ref_in = pos[sel][:, out_c].mean(axis=1) - pos[sel][:, in_c].mean(axis=1)
            if len(in_c) == 1 or len(out_c) == 1:
                lone = in_c if len(in_c) == 1 else out_c
                rest = out_c if len(in_c) == 1 else in_c
                e = [edge(lone[0], r, sel) for r in rest]
                tris_keys.append(np.stack([k for k, _ in e], axis=1))
                edge_key += [k for k, _ in e]
                edge_pts += [p for _, p in e]
                tris_inside.append(ref_in)
            else:
                a, b = in_c
                c1, d1 = out_c
                e = [edge(a, c1, sel), edge(a, d1, sel), edge(b, d1, sel), edge(b, c1, sel)]
                k = [kk_ for kk_, _ in e]
                edge_key += k
                edge_pts += [p for _, p in e]
                tris_keys.append(np.stack([k[0], k[1], k[2]], axis=1))
                tris_keys.append(np.stack([k[0], k[2], k[3]], axis=1))
                tris_inside += [ref_in, ref_in]

    if not tris_keys:
        return np.zeros((0, 3)), np.zeros((0, 3), int), np.zeros(0, np.int64), nx * ny * nz
    all_keys = np.concatenate(edge_key)
    all_pts = np.concatenate(edge_pts)
    uniq, first = np.unique(all_keys, return_index=True)
    verts = all_pts[first]
    tri_keys = np.concatenate(tris_keys)
    tris = np.searchsorted(uniq, tri_keys)
    ref = np.concatenate(tris_inside)

    # drop triangles collapsed by snapping, and duplicate (coincident) faces
    good = (tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 0] != tris[:, 2])
    tris, ref = tris[good], ref[good]
    _, inv, cnt = np.unique(np.sort(tris, axis=1), axis=0, return_inverse=True, return_counts=True)
    keep = cnt[inv.ravel()] == 1  # coincident faces from snapping cancel each other
    tris, ref = tris[keep], ref[keep]
    used, tris = np.unique(tris, return_inverse=True)
    tris = tris.reshape(-1, 3)
    verts, uniq = verts[used], uniq[used]
    # orient normals outward (away from the inside corner)
    v0, v1, v2 = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    nrm = np.cross(v1 - v0, v2 - v0)
    flip = (nrm * ref).sum(axis=1) < 0
    tris[flip] = tris[flip][:, ::-1]
    verts = origin + verts * spacing
    return verts, tris, uniq, nx * ny * nz


def mesh_volume(verts, tris):
    v0, v1, v2 = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    c = verts.mean(axis=0)
    return float(np.einsum("ij,ij->i", v0 - c, np.cross(v1 - c, v2 - c)).sum() / 6.0)


def mesh_area(verts, tris):
    v0, v1, v2 = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    return float(np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1).sum() / 2.0)


def is_closed(tris):
    e = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return bool((counts == 2).all())


def grid_surface(xs, ys, zgrid):
    """Triangulate a regular-grid DTM; zgrid has shape (len(xs), len(ys))."""
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    verts = np.stack([X.ravel(), Y.ravel(), zgrid.ravel()], axis=1)
    ny = len(ys)
    i, j = np.meshgrid(np.arange(len(xs) - 1), np.arange(ny - 1), indexing="ij")
    a = (i * ny + j).ravel()
    b, c, d = a + ny, a + ny + 1, a + 1
    tris = np.concatenate([np.stack([a, b, c], 1), np.stack([a, c, d], 1)])
    return verts, tris


def marching_squares(f, iso=0.0):
    """Contour a 2D field (positive inside); returns list of (points[n,2] in index
    space, closed flag).  The field is padded so contours are always closed."""
    f = np.pad(np.asarray(f, float) - iso, 1, constant_values=-1.0)
    f[f == 0.0] = 1e-9
    nx, ny = f.shape
    ins = f > 0

    def ekey(i0, j0, i1, j1):
        return (min(i0 * ny + j0, i1 * ny + j1), max(i0 * ny + j0, i1 * ny + j1))

    def ept(i0, j0, i1, j1):
        a, b = f[i0, j0], f[i1, j1]
        t = a / (a - b)
        return (i0 + t * (i1 - i0) - 1.0, j0 + t * (j1 - j0) - 1.0)

    pts, adj = {}, {}

    def link(k1, k2):
        adj.setdefault(k1, []).append(k2)
        adj.setdefault(k2, []).append(k1)

    I, J = np.nonzero(ins[:-1, :-1] ^ ins[1:, :-1] | ins[:-1, :-1] ^ ins[:-1, 1:] | ins[1:, 1:] ^ ins[1:, :-1])
    for i, j in zip(I.tolist(), J.tolist()):
        c = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        e = []
        for n in range(4):
            p, q = c[n], c[(n + 1) % 4]
            if ins[p] != ins[q]:
                k = ekey(*p, *q)
                pts[k] = ept(*p, *q)
                e.append((n, k))
        if len(e) == 2:
            link(e[0][1], e[1][1])
        elif len(e) == 4:  # saddle: resolve with the cell-centre average
            centre_in = f[i, j] + f[i + 1, j] + f[i + 1, j + 1] + f[i, j + 1] > 0
            if centre_in == bool(ins[c[0]]):
                link(e[0][1], e[1][1]); link(e[2][1], e[3][1])
            else:
                link(e[0][1], e[3][1]); link(e[1][1], e[2][1])

    lines, seen = [], set()
    for start in adj:
        if start in seen:
            continue
        chain, prev, cur = [start], None, start
        seen.add(start)
        while True:
            nxt = [k for k in adj[cur] if k != prev and k not in seen]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            seen.add(cur)
            chain.append(cur)
        closed = start in adj[cur] and len(chain) > 2
        lines.append((np.array([pts[k] for k in chain]), closed))
    return lines
