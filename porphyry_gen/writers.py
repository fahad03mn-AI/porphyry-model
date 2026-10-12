"""Writers for common mining-software exchange formats (OBJ, DXF R12, Surpac STR, CSV)."""
import csv
import datetime as _dt


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def write_obj(path, verts, tris, name):
    with open(path, "w") as fh:
        fh.write(f"# {name}\n# vertices: {len(verts)}  triangles: {len(tris)}\no {name}\n")
        fh.writelines(f"v {x:.2f} {y:.2f} {z:.2f}\n" for x, y, z in verts)
        fh.writelines(f"f {a + 1} {b + 1} {c + 1}\n" for a, b, c in tris)


def _dxf_open(fh):
    fh.write("0\nSECTION\n2\nENTITIES\n")


def _dxf_close(fh):
    fh.write("0\nENDSEC\n0\nEOF\n")


def write_dxf_mesh(path, verts, tris, layer):
    """Triangles as DXF R12 3DFACE entities (imports into Surpac, Datamine,
    Vulcan, Micromine, Leapfrog, Deswik, AutoCAD...)."""
    with open(path, "w") as fh:
        _dxf_open(fh)
        for t in tris:
            p = verts[[t[0], t[1], t[2], t[2]]]
            fh.write(f"0\n3DFACE\n8\n{layer}\n")
            for n, (x, y, z) in enumerate(p):
                fh.write(f"1{n}\n{x:.2f}\n2{n}\n{y:.2f}\n3{n}\n{z:.2f}\n")
        _dxf_close(fh)


def write_dxf_polylines(path, strings):
    """strings: iterable of (layer, points[n,3], closed)."""
    with open(path, "w") as fh:
        _dxf_open(fh)
        for layer, pts, closed in strings:
            fh.write(f"0\nPOLYLINE\n8\n{layer}\n66\n1\n10\n0.0\n20\n0.0\n30\n0.0\n70\n{9 if closed else 8}\n")
            for x, y, z in pts:
                fh.write(f"0\nVERTEX\n8\n{layer}\n10\n{x:.3f}\n20\n{y:.3f}\n30\n{z:.3f}\n70\n32\n")
            fh.write(f"0\nSEQEND\n8\n{layer}\n")
        _dxf_close(fh)


def write_surpac_str(path, title, strings):
    """strings: iterable of (string_no, points[n,3], closed, description).
    Surpac STR stores Y (northing) before X (easting)."""
    date = _dt.date(2026, 1, 1).strftime("%d-%b-%y")
    with open(path, "w") as fh:
        fh.write(f"{title},{date},,ssi_styles:styles.ssi\n")
        fh.write("0, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000\n")
        for sno, pts, closed, desc in strings:
            seq = list(pts) + ([pts[0]] if closed else [])
            for x, y, z in seq:
                fh.write(f"{sno}, {y:.3f}, {x:.3f}, {z:.3f}, {desc}\n")
            fh.write("0, 0.000, 0.000, 0.000,\n")
        fh.write("0, 0.000, 0.000, 0.000, END\n")
