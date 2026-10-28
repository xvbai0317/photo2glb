# -*- coding: utf-8 -*-
"""
Procedural rebuild of the Parker Hannifin 3-valve instrument manifold
(angled-pattern block valve) from the supplied reference photo.

Authored in millimetres, Y up, +Z toward the viewer (front face), +X right;
exported as a metric GLB (0.001 scale).

Run:  python build_manifold.py
"""
import json
import math
import os
import struct

import numpy as np
from PIL import Image

TAU = 2.0 * math.pi
D2R = math.pi / 180.0
WORK = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = WORK
REF = os.path.join(os.path.dirname(WORK), "docs", "reference.png")


# ---------------------------------------------------------------- math utils
def T(x, y, z):
    m = np.eye(4)
    m[:3, 3] = (x, y, z)
    return m


def Rx(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, 1]], float)


def Rz(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0, 0], [s, c, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], float)


def normalize(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def inset_poly(poly, d):
    """Miter inset of a convex CCW polygon by distance d."""
    p = np.asarray(poly, float)
    n = len(p)
    out = []
    for i in range(n):
        a, b, c = p[(i - 1) % n], p[i], p[(i + 1) % n]
        e0 = normalize(b - a)
        e1 = normalize(c - b)
        in0 = np.array([-e0[1], e0[0]])
        in1 = np.array([-e1[1], e1[0]])
        A = np.array([[e0[0], -e1[0]], [e0[1], -e1[1]]])
        rhs = (b + in1 * d) - (b + in0 * d)
        try:
            t = np.linalg.solve(A, rhs)[0]
        except np.linalg.LinAlgError:
            t = 0.0
        out.append(b + in0 * d + e0 * t)
    return np.array(out)


def ray_poly_hit(poly, origin, ang):
    d = np.array([math.cos(ang), math.sin(ang)])
    p = np.asarray(poly, float)
    best = None
    for i in range(len(p)):
        a, b = p[i], p[(i + 1) % len(p)]
        e = b - a
        den = d[0] * (-e[1]) - d[1] * (-e[0])
        if abs(den) < 1e-12:
            continue
        w = a - origin
        t = (w[0] * (-e[1]) - w[1] * (-e[0])) / den
        s = (d[0] * w[1] - d[1] * w[0]) / den
        if t > 1e-9 and -1e-9 <= s <= 1.0 + 1e-9:
            if best is None or t < best:
                best = t
    return best if best is not None else 1.0


# ---------------------------------------------------------------- geo buffer
class Geo:
    """Accumulates triangles with per-vertex normals and optional UVs."""

    def __init__(self):
        self.P = []
        self.N = []
        self.UV = []
        self.T = []
        self.has_uv = False

    def flat_quad(self, pts, outward, uv=None):
        p = np.asarray(pts, float)
        n = np.cross(p[1] - p[0], p[2] - p[0])
        if np.linalg.norm(n) < 1e-9:
            n = np.cross(p[2] - p[0], p[3] - p[0])
        n = normalize(n)
        order = [0, 1, 2, 3]
        if np.dot(n, np.asarray(outward, float)) < 0:
            order = [0, 3, 2, 1]
            n = -n
        base = len(self.P)
        for k in order:
            self.P.append(p[k])
            self.N.append(n)
            self.UV.append(uv[k] if uv is not None else (0.0, 0.0))
        self.T += [(base, base + 1, base + 2), (base, base + 2, base + 3)]
        if uv is not None:
            self.has_uv = True

    def flat_tri(self, pts, outward, uv=None):
        p = np.asarray(pts, float)
        n = normalize(np.cross(p[1] - p[0], p[2] - p[0]))
        order = [0, 1, 2]
        if np.dot(n, np.asarray(outward, float)) < 0:
            order = [0, 2, 1]
            n = -n
        base = len(self.P)
        for k in order:
            self.P.append(p[k])
            self.N.append(n)
            self.UV.append(uv[k] if uv is not None else (0.0, 0.0))
        self.T.append((base, base + 1, base + 2))
        if uv is not None:
            self.has_uv = True

    def smooth_quad(self, pts, nrms, outward, uv=None):
        p = np.asarray(pts, float)
        nn = np.asarray(nrms, float)
        order = [0, 1, 2, 3]
        ngeo = np.cross(p[1] - p[0], p[2] - p[0])
        if np.linalg.norm(ngeo) < 1e-12:
            ngeo = np.cross(p[2] - p[0], p[3] - p[0])
        # keep the winding consistent with the stored (outward) vertex normals
        if np.dot(ngeo, nn.sum(axis=0)) < 0:
            order = [0, 3, 2, 1]
        base = len(self.P)
        for k in order:
            self.P.append(p[k])
            self.N.append(normalize(nn[k]))
            self.UV.append(uv[k] if uv is not None else (0.0, 0.0))
        self.T += [(base, base + 1, base + 2), (base, base + 2, base + 3)]
        if uv is not None:
            self.has_uv = True

    def extend(self, other, matrix=None):
        base = len(self.P)
        if matrix is None:
            self.P += list(other.P)
            self.N += list(other.N)
        else:
            R, t = matrix[:3, :3], matrix[:3, 3]
            Nm = np.linalg.inv(R).T
            for p in other.P:
                self.P.append(R.dot(p) + t)
            for n in other.N:
                self.N.append(normalize(Nm.dot(n)))
        self.UV += list(other.UV)
        self.T += [(a + base, b + base, c + base) for (a, b, c) in other.T]
        self.has_uv = self.has_uv or other.has_uv
        return self

    def arrays(self):
        if not self.P:
            return (np.zeros((0, 3)), np.zeros((0, 3)), np.zeros((0, 2)), np.zeros((0, 3), np.int64))
        return (np.array(self.P, float).reshape(-1, 3),
                np.array(self.N, float).reshape(-1, 3),
                np.array(self.UV, float).reshape(-1, 2),
                np.array(self.T, np.int64).reshape(-1, 3))


# ---------------------------------------------------------------- generators
def g_revolve(profile, seg=64, smooth=True, cap_top=False, cap_bottom=False, flip=False):
    g = Geo()
    n = len(profile)
    segn = []
    for i in range(n - 1):
        r0, y0 = profile[i]
        r1, y1 = profile[i + 1]
        dr, dy = r1 - r0, y1 - y0
        L = math.hypot(dr, dy)
        segn.append((dy / L, -dr / L))
    if smooth:
        pn = []
        for i in range(n):
            if i == 0:
                pn.append(segn[0])
            elif i == n - 1:
                pn.append(segn[-1])
            else:
                pn.append((segn[i - 1][0] + segn[i][0], segn[i - 1][1] + segn[i][1]))
        rings = []
        for i in range(n):
            r, y = profile[i]
            nr, ny = normalize(pn[i])
            ring = []
            for j in range(seg):
                a = j / seg * TAU
                ring.append(((r * math.cos(a), y, r * math.sin(a)),
                             (nr * math.cos(a), ny, nr * math.sin(a))))
            rings.append(ring)
        for i in range(n - 1):
            for j in range(seg):
                j2 = (j + 1) % seg
                p = [rings[i][j][0], rings[i][j2][0], rings[i + 1][j2][0], rings[i + 1][j][0]]
                nn = [rings[i][j][1], rings[i][j2][1], rings[i + 1][j2][1], rings[i + 1][j][1]]
                g.smooth_quad(p, nn, (p[0][0], 0, p[0][2]))
    else:
        for i in range(n - 1):
            r0, y0 = profile[i]
            r1, y1 = profile[i + 1]
            nr, ny = segn[i]
            for j in range(seg):
                a0 = j / seg * TAU
                a1 = (j + 1) / seg * TAU
                p = [(r0 * math.cos(a0), y0, r0 * math.sin(a0)),
                     (r0 * math.cos(a1), y0, r0 * math.sin(a1)),
                     (r1 * math.cos(a1), y1, r1 * math.sin(a1)),
                     (r1 * math.cos(a0), y1, r1 * math.sin(a0))]
                nn = [(nr * math.cos(a0), ny, nr * math.sin(a0)),
                      (nr * math.cos(a1), ny, nr * math.sin(a1)),
                      (nr * math.cos(a1), ny, nr * math.sin(a1)),
                      (nr * math.cos(a0), ny, nr * math.sin(a0))]
                g.smooth_quad(p, nn, (math.cos(a0), 0, math.sin(a0)))
    if cap_top:
        r, y = profile[-1]
        for j in range(seg):
            a0 = j / seg * TAU
            a1 = (j + 1) / seg * TAU
            g.flat_tri([(0, y, 0), (r * math.cos(a0), y, r * math.sin(a0)), (r * math.cos(a1), y, r * math.sin(a1))],
                       (0, 1, 0))
    if cap_bottom:
        r, y = profile[0]
        for j in range(seg):
            a0 = j / seg * TAU
            a1 = (j + 1) / seg * TAU
            g.flat_tri([(0, y, 0), (r * math.cos(a0), y, r * math.sin(a0)), (r * math.cos(a1), y, r * math.sin(a1))],
                       (0, -1, 0))
    if flip:
        g.N = [-np.asarray(x) for x in g.N]
        g.T = [(a, c, b) for (a, b, c) in g.T]
    return g


def g_cyl(r, y0, y1, seg=64, chamfer_top=0.0, chamfer_bot=0.0, cap=True):
    prof = []
    if chamfer_bot > 0:
        prof += [(r - chamfer_bot, y0), (r, y0 + chamfer_bot)]
    else:
        prof += [(r, y0)]
    if chamfer_top > 0:
        prof += [(r, y1 - chamfer_top), (r - chamfer_top, y1)]
    else:
        prof += [(r, y1)]
    return g_revolve(prof, seg=seg, smooth=True, cap_top=cap, cap_bottom=cap)


def g_ngon(n, af, y0, y1, phase=0.0, chamfer_top=0.0, chamfer_bot=0.0, cap=True):
    r = (af / 2.0) / math.cos(math.pi / n)
    rings = []
    if chamfer_bot > 0:
        rings += [(r - chamfer_bot, y0), (r, y0 + chamfer_bot)]
    else:
        rings.append((r, y0))
    if chamfer_top > 0:
        rings += [(r, y1 - chamfer_top), (r - chamfer_top, y1)]
    else:
        rings.append((r, y1))
    g = Geo()
    angs = [phase + k * TAU / n for k in range(n)]
    for ri in range(len(rings) - 1):
        r0, y0_ = rings[ri]
        r1, y1_ = rings[ri + 1]
        for k in range(n):
            a0 = angs[k] - math.pi / n
            a1 = angs[k] + math.pi / n
            p = [(r0 * math.cos(a0), y0_, r0 * math.sin(a0)),
                 (r0 * math.cos(a1), y0_, r0 * math.sin(a1)),
                 (r1 * math.cos(a1), y1_, r1 * math.sin(a1)),
                 (r1 * math.cos(a0), y1_, r1 * math.sin(a0))]
            g.flat_quad(p, (math.cos(angs[k]), 0, math.sin(angs[k])))
    if cap:
        for (rr, yy, dirn) in ((rings[-1][0], rings[-1][1], 1), (rings[0][0], rings[0][1], -1)):
            for k in range(n):
                a0 = angs[k] - math.pi / n
                a1 = angs[k] + math.pi / n
                g.flat_tri([(0, yy, 0), (rr * math.cos(a0), yy, rr * math.sin(a0)),
                            (rr * math.cos(a1), yy, rr * math.sin(a1))], (0, dirn, 0))
    return g


def g_fluted(r, y0, y1, flutes=18, amp=0.4, seg_per=7, taper=0.55):
    g = Geo()
    seg = flutes * seg_per
    levels = [(y0, r - taper), (y0 + 0.5, r), (y1 - 0.5, r), (y1, r - taper)]
    rings = []
    for (y, rr) in levels:
        ring = []
        for j in range(seg):
            a = j / seg * TAU
            rad = rr + amp * math.cos(flutes * a)
            ring.append((rad * math.cos(a), y, rad * math.sin(a)))
        rings.append(ring)
    for i in range(len(rings) - 1):
        for j in range(seg):
            j2 = (j + 1) % seg
            a0 = j / seg * TAU
            g.flat_quad([rings[i][j], rings[i][j2], rings[i + 1][j2], rings[i + 1][j]],
                        (math.cos(a0), 0.0, math.sin(a0)))
    for (ring, y, dirn) in ((rings[-1], y1, 1), (rings[0], y0, -1)):
        for j in range(seg):
            j2 = (j + 1) % seg
            g.flat_tri([(0, y, 0), ring[j], ring[j2]], (0, dirn, 0))
    return g


def g_threaded(r, y0, y1, rings=5, amp=0.35, taper=0.0, seg=56, cap_top=True, cap_bottom=True):
    prof = [(r, y0)]
    h = y1 - y0
    steps = rings * 4
    for k in range(1, steps + 1):
        t = k / steps
        y = y0 + h * t
        rr = r - taper * t + amp * (0.5 + 0.5 * math.cos(TAU * rings * t)) ** 2
        prof.append((rr, y))
    return g_revolve(prof, seg=seg, smooth=False, cap_top=cap_top, cap_bottom=cap_bottom)


def g_bore(r_mouth, r_wall, depth, rings=6, amp=0.45, taper=0.5, seg=64):
    prof = [(r_mouth, 0.0), (r_wall, 1.3)]
    steps = rings * 4
    for k in range(1, steps + 1):
        t = k / steps
        y = 1.3 + (depth - 1.3) * t
        rr = r_wall - taper * t + amp * (0.5 + 0.5 * math.cos(TAU * rings * t)) ** 2
        prof.append((rr, y))
    g = g_revolve(prof, seg=seg, smooth=False, cap_bottom=False, cap_top=False, flip=True)
    r_end = prof[-1][0]
    for j in range(seg):
        a0 = j / seg * TAU
        a1 = (j + 1) / seg * TAU
        g.flat_tri([(0, depth, 0), (r_end * math.cos(a0), depth, r_end * math.sin(a0)),
                    (r_end * math.cos(a1), depth, r_end * math.sin(a1))], (0, -1, 0))
    return g


def emit_face_hole(g, origin, u, v, quad2d, hole_c, hole_r, outward, seg=72):
    o, u, v = np.asarray(origin, float), np.asarray(u, float), np.asarray(v, float)
    out = np.asarray(outward, float)
    hc = np.asarray(hole_c, float)

    def p3(pt):
        return o + u * pt[0] + v * pt[1]

    for j in range(seg):
        a0, a1 = j / seg * TAU, (j + 1) / seg * TAU
        i0 = hc + hole_r * np.array([math.cos(a0), math.sin(a0)])
        i1 = hc + hole_r * np.array([math.cos(a1), math.sin(a1)])
        o0 = hc + ray_poly_hit(quad2d, hc, a0) * np.array([math.cos(a0), math.sin(a0)])
        o1 = hc + ray_poly_hit(quad2d, hc, a1) * np.array([math.cos(a1), math.sin(a1)])
        g.flat_quad([p3(i0), p3(i1), p3(o1), p3(o0)], out)


def g_wheel(L=33.0, W=16.0, thick=6.4, hub_r=5.0,
            hole_dx=11.3, hole_r=3.05, stem_r=2.9,
            lobes=6, amp=0.085, nx=76, bev=0.85):
    """Parker-style flat handwheel: an elongated scalloped plate with two round
    finger holes plus a centre hole for the stem, bevelled rim, flat shaded."""
    R = W / 2.0
    st = L / 2.0 - R
    holes = [(-hole_dx, hole_r), (0.0, stem_r), (hole_dx, hole_r)]

    def zout(x):
        ax = abs(x)
        if ax <= st:
            z = R
        else:
            u = min(ax, L / 2.0) - st
            z = math.sqrt(max(0.0, R * R - u * u))
        lim = 0.60 * L / 2.0
        fade = 1.0
        if ax > lim:
            t = min(1.0, (ax - lim) / (0.40 * L / 2.0))
            fade = 1.0 - t * t * (3.0 - 2.0 * t)
        z *= 1.0 + amp * fade * math.cos(TAU * (lobes / 2.0) * x / L)
        return max(z, 0.05)

    def h_at(x):
        m = 0.0
        for (x0, r) in holes:
            d = r * r - (x - x0) ** 2
            if d > 0:
                m = max(m, math.sqrt(d))
        return m

    xs = [-L / 2.0 + i * L / nx for i in range(nx + 1)]
    ytop, ybot = thick / 2.0, -thick / 2.0
    g = Geo()
    for i in range(nx):
        xa, xb = xs[i], xs[i + 1]
        za, zb = zout(xa), zout(xb)
        ha, hb = h_at(xa), h_at(xb)
        if ha <= 1e-9 and hb <= 1e-9:
            g.flat_quad([(xa, ytop, za), (xb, ytop, zb), (xb, ytop, -zb), (xa, ytop, -za)], (0, 1, 0))
            g.flat_quad([(xa, ybot, -za), (xb, ybot, -zb), (xb, ybot, zb), (xa, ybot, za)], (0, -1, 0))
        else:
            g.flat_quad([(xa, ytop, za), (xb, ytop, zb), (xb, ytop, hb), (xa, ytop, ha)], (0, 1, 0))
            g.flat_quad([(xa, ytop, -ha), (xb, ytop, -hb), (xb, ytop, -zb), (xa, ytop, -za)], (0, 1, 0))
            g.flat_quad([(xa, ybot, ha), (xb, ybot, hb), (xb, ybot, zb), (xa, ybot, za)], (0, -1, 0))
            g.flat_quad([(xa, ybot, -za), (xb, ybot, -zb), (xb, ybot, -hb), (xa, ybot, -ha)], (0, -1, 0))
            xm = (xa + xb) / 2.0
            xc = min((h[0] for h in holes), key=lambda c: abs(c - xm))
            nrm = normalize(np.array([-(hb - ha), 0.0, xb - xa]))
            if np.dot(nrm, np.array([xm - xc, 0.0, (ha + hb) / 2.0])) > 0:
                nrm = -nrm
            g.flat_quad([(xa, ytop, ha), (xb, ytop, hb), (xb, ybot, hb), (xa, ybot, ha)], nrm)
            g.flat_quad([(xa, ytop, -ha), (xb, ytop, -hb), (xb, ybot, -hb), (xa, ybot, -ha)],
                        (nrm[0], 0.0, -nrm[2]))
        for sgn in (1.0, -1.0):
            A0 = np.array([xa, za * sgn])
            A1 = np.array([xb, zb * sgn])
            e = A1 - A0
            n2 = normalize(np.array([-e[1], e[0]]))
            if np.dot(n2, (A0 + A1)) < 0:
                n2 = -n2
            out = np.array([n2[0], 0.0, n2[1]])
            O0 = np.array([xa, 0.0, za * sgn])
            O1 = np.array([xb, 0.0, zb * sgn])
            I0 = O0 - out * bev
            I1 = O1 - out * bev
            g.flat_quad([O0 + [0, ytop, 0], O1 + [0, ytop, 0],
                         I1 + [0, ytop - bev, 0], I0 + [0, ytop - bev, 0]], out)
            g.flat_quad([I0 + [0, ytop - bev, 0], I1 + [0, ytop - bev, 0],
                         I1 + [0, ybot + bev, 0], I0 + [0, ybot + bev, 0]], out)
            g.flat_quad([I0 + [0, ybot + bev, 0], I1 + [0, ybot + bev, 0],
                         O1 + [0, ybot, 0], O0 + [0, ybot, 0]], out)

    g.extend(g_cyl(hub_r, ytop, ytop + 3.2, seg=56, chamfer_top=0.9))
    g.extend(g_cyl(hub_r, ybot - 1.5, ybot, seg=56, chamfer_bot=0.7))
    return g


# ---------------------------------------------------------------- materials
def lin(c):
    return [round(max(0.0, min(1.0, v)) ** 2.2, 5) for v in c]


MATERIALS = [
    {"name": "StainlessSteel", "pbr": {"baseColorFactor": lin((0.90, 0.905, 0.915)) + [1.0], "metallicFactor": 0.95, "roughnessFactor": 0.30}},
    {"name": "StainlessMachined", "pbr": {"baseColorFactor": lin((0.89, 0.895, 0.905)) + [1.0], "metallicFactor": 0.95, "roughnessFactor": 0.36}},
    {"name": "ThreadCavity", "pbr": {"baseColorFactor": lin((0.30, 0.31, 0.33)) + [1.0], "metallicFactor": 0.7, "roughnessFactor": 0.50}},
    {"name": "RedPolymer", "pbr": {"baseColorFactor": lin((0.86, 0.16, 0.11)) + [1.0], "metallicFactor": 0.0, "roughnessFactor": 0.42}},
    {"name": "BluePolymer", "pbr": {"baseColorFactor": lin((0.10, 0.27, 0.72)) + [1.0], "metallicFactor": 0.0, "roughnessFactor": 0.40}},
    {"name": "NamePlate", "pbr": {"baseColorFactor": [1.0, 1.0, 1.0, 1.0], "metallicFactor": 0.0, "roughnessFactor": 0.55}},
]
M_STEEL, M_MACH, M_CAV, M_RED, M_BLUE, M_PLATE = range(6)

POLY = [(-33.0, -20.0), (33.0, -20.0), (33.0, 6.5), (13.0, 20.0), (-13.0, 20.0), (-33.0, 6.5)]
Z_F, Z_B, BEV = 20.0, -20.0, 1.6


# ---------------------------------------------------------------- label texture
def make_label_texture():
    src = Image.open(REF).convert("RGB")
    quad = (345.0, 513.0, 348.5, 647.0, 515.0, 655.0, 511.5, 523.0)  # TL, BL, BR, TR
    txt = src.transform((420, 360), Image.QUAD, data=quad, resample=Image.BICUBIC)
    a = np.asarray(txt, float) / 255.0
    lum = a.mean(axis=2)
    white = float(np.percentile(lum, 88))
    ink = float(np.percentile(lum, 2))
    v = np.clip((lum - ink * 0.55) / max(1e-3, white - ink * 0.55), 0, 1)
    v = 0.06 + 0.90 * v ** 1.05
    rgb = np.stack([v * 0.99, v * 0.99, v * 0.975], axis=2)
    img = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8))
    path = os.path.join(WORK, "label.png")
    img.save(path)
    return path, img.size


# ---------------------------------------------------------------- assembly
class Node:
    def __init__(self, name, matrix=None):
        self.name = name
        self.local = np.eye(4) if matrix is None else matrix
        self.children = []
        self.parts = []   # (material_index, Geo)

    def add(self, mat, geo, matrix=None):
        if matrix is not None:
            geo = Geo().extend(geo, matrix)
        if len(geo.P):
            self.parts.append((mat, geo))
        return self

    def world(self, parent=np.eye(4)):
        w = parent.dot(self.local)
        return w


def build_body(node):
    ins = inset_poly(POLY, BEV)
    g_front = Geo()
    emit_face_hole(g_front, (0, 0, Z_F), (1, 0, 0), (0, 1, 0), ins, (-1.5, -8.0), 5.8, (0, 0, 1), seg=72)
    node.add(M_MACH, g_front)

    g_back = Geo()
    c = np.mean(ins, axis=0)
    for i in range(len(ins)):
        a, b = ins[i], ins[(i + 1) % len(ins)]
        g_back.flat_tri([(c[0], c[1], Z_B), (a[0], a[1], Z_B), (b[0], b[1], Z_B)], (0, 0, -1))
    node.add(M_MACH, g_back)

    for i in range(len(POLY)):
        a, b = np.array(POLY[i]), np.array(POLY[(i + 1) % len(POLY)])
        e = b - a
        outn = normalize(np.array([e[1], -e[0], 0.0]))
        g = Geo()
        g.flat_quad([(a[0], a[1], Z_F - BEV), (b[0], b[1], Z_F - BEV),
                     (ins[(i + 1) % len(ins)][0], ins[(i + 1) % len(ins)][1], Z_F),
                     (ins[i][0], ins[i][1], Z_F)], outn + np.array([0, 0, 1.0]))
        g.flat_quad([(b[0], b[1], Z_B + BEV), (a[0], a[1], Z_B + BEV),
                     (ins[i][0], ins[i][1], Z_B),
                     (ins[(i + 1) % len(ins)][0], ins[(i + 1) % len(ins)][1], Z_B)], outn + np.array([0, 0, -1.0]))
        node.add(M_MACH, g)

        gw = Geo()
        if i == 1:
            emit_face_hole(gw, (33.0, -20.0, -18.4), (0, 1, 0), (0, 0, 1),
                           [(0, 0), (26.5, 0), (26.5, 36.8), (0, 36.8)], (13.0, 18.4), 7.0, (1, 0, 0), seg=72)
        elif i == 5:
            emit_face_hole(gw, (-33.0, -20.0, -18.4), (0, 1, 0), (0, 0, 1),
                           [(0, 0), (26.5, 0), (26.5, 36.8), (0, 36.8)], (13.0, 18.4), 7.0, (-1, 0, 0), seg=72)
        else:
            gw.flat_quad([(a[0], a[1], Z_B + BEV), (b[0], b[1], Z_B + BEV),
                          (b[0], b[1], Z_F - BEV), (a[0], a[1], Z_F - BEV)], outn)
        node.add(M_MACH, gw)

        if i == 1:
            node.add(M_CAV, g_bore(7.0, 5.4, 14.0, rings=6, amp=0.4, taper=0.4, seg=64),
                     T(33.0, -7.0, 0.0).dot(Rz(90 * D2R)))
        elif i == 5:
            node.add(M_CAV, g_bore(7.0, 5.4, 14.0, rings=6, amp=0.4, taper=0.4, seg=64),
                     T(-33.0, -7.0, 0.0).dot(Rz(-90 * D2R)))

    node.add(M_CAV, g_bore(5.8, 4.4, 12.0, rings=6, amp=0.35, taper=0.35, seg=64),
             T(-1.5, -8.0, Z_F).dot(Rx(-90 * D2R)))


def build_label(node):
    w, h = 20.0, 17.0
    cx, cy, ang = -6.0, 9.0, 1.2 * D2R
    pos, uv = [], [(0, 1), (1, 1), (1, 0), (0, 0)]
    for (sx, sy) in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        x, y = sx * w / 2, sy * h / 2
        xr = x * math.cos(ang) - y * math.sin(ang)
        yr = x * math.sin(ang) + y * math.cos(ang)
        pos.append((cx + xr, cy + yr, Z_F + 0.03))
    g = Geo()
    g.flat_quad(pos, (0, 0, 1), uv)
    node.add(M_PLATE, g)


def build_valve(node, s, band_mat, wheel_scale):
    g = Geo()
    g.extend(g_ngon(6, 20.0 * s, 0.0, 5.0 * s, chamfer_top=1.1 * s))
    g.extend(g_cyl(7.5 * s, 4.0 * s, 11.6 * s, seg=56, chamfer_top=1.0 * s, chamfer_bot=0.6 * s))
    g.extend(g_ngon(6, 15.0 * s, 11.4 * s, 15.4 * s, chamfer_top=0.7 * s))
    g.extend(g_threaded(5.3 * s, 15.2 * s, 18.4 * s, rings=4, amp=0.22 * s, seg=48))
    g.extend(g_ngon(6, 11.6 * s, 18.2 * s, 22.0 * s, chamfer_top=0.6 * s))
    node.add(M_STEEL, g)

    gb = Geo()
    gb.extend(g_fluted(6.9 * s, 21.9 * s, 25.4 * s, flutes=16, amp=0.35 * s, seg_per=6, taper=0.5 * s))
    gb.extend(g_cyl(7.15 * s, 25.2 * s, 26.0 * s, seg=48, chamfer_top=0.4 * s))
    node.add(band_mat, gb)

    gtop = Geo()
    gtop.extend(g_ngon(6, 12.2 * s, 25.8 * s, 29.0 * s, chamfer_top=0.7 * s))
    gtop.extend(g_cyl(2.15 * s, 28.8 * s, 34.2 * s, seg=24))
    node.add(M_STEEL, gtop)

    ws = wheel_scale
    wy = 34.4 * s
    wheel_geo = Geo().extend(g_wheel(), T(0, wy, 0).dot(np.diag([ws, ws, ws, 1.0])))
    node.add(M_STEEL, wheel_geo)


# ---------------------------------------------------------------- GLB writer
def write_glb(path, root, materials, texture_png=None, root_matrix=None):
    order = []

    def walk(n):
        order.append(n)
        for c in n.children:
            walk(c)
    walk(root)

    bin_parts, buffer_views, accessors, meshes = [], [], [], []
    offset = [0]
    img_view = [None]

    def add_view(data, target=None):
        while offset[0] % 4:
            bin_parts.append(b"\x00")
            offset[0] += 1
        buffer_views.append({"buffer": 0, "byteOffset": offset[0], "byteLength": len(data),
                             **({"target": target} if target else {})})
        bin_parts.append(data)
        offset[0] += len(data)
        return len(buffer_views) - 1

    if texture_png is not None:
        img_view[0] = add_view(open(texture_png, "rb").read())

    nodes_json = []
    index = {}
    for n in order:
        prims = []
        for (mat, geo) in n.parts:
            P, Nn, UV, Tr = geo.arrays()
            if len(P) == 0:
                continue
            attrs = {}
            pv = add_view(P.astype(np.float32).tobytes(), 34962)
            accessors.append({"bufferView": pv, "componentType": 5126, "count": int(len(P)),
                              "type": "VEC3",
                              "min": [round(float(x), 5) for x in P.min(axis=0)],
                              "max": [round(float(x), 5) for x in P.max(axis=0)]})
            attrs["POSITION"] = len(accessors) - 1
            nv = add_view(Nn.astype(np.float32).tobytes(), 34962)
            accessors.append({"bufferView": nv, "componentType": 5126, "count": int(len(Nn)), "type": "VEC3"})
            attrs["NORMAL"] = len(accessors) - 1
            if geo.has_uv:
                tv = add_view(UV.astype(np.float32).tobytes(), 34962)
                accessors.append({"bufferView": tv, "componentType": 5126, "count": int(len(UV)), "type": "VEC2"})
                attrs["TEXCOORD_0"] = len(accessors) - 1
            iv = add_view(Tr.astype(np.uint32).tobytes(), 34963)
            accessors.append({"bufferView": iv, "componentType": 5125, "count": int(Tr.size), "type": "SCALAR"})
            prims.append({"attributes": attrs, "indices": len(accessors) - 1, "material": mat, "mode": 4})
        meshes.append({"primitives": prims, "name": n.name}) if prims else None
        lm = root_matrix if (root_matrix is not None and n is root) else n.local
        entry = {"name": n.name, "matrix": lm.T.flatten().round(9).tolist()}
        if prims:
            entry["mesh"] = len(meshes) - 1
        index[id(n)] = len(nodes_json)
        nodes_json.append(entry)
    for n in order:
        if n.children:
            nodes_json[index[id(n)]]["children"] = [index[id(c)] for c in n.children]

    while offset[0] % 4:
        bin_parts.append(b"\x00")
        offset[0] += 1
    BIN = b"".join(bin_parts)

    mats_json = []
    for m in materials:
        mats_json.append({"name": m["name"], "pbrMetallicRoughness": m["pbr"],
                          "doubleSided": bool(m.get("doubleSided", False))})

    gltf = {
        "asset": {"version": "2.0", "generator": "img2threejs procedural manifold rebuild"},
        "scene": 0,
        "scenes": [{"name": "ParkerThreeValveManifold", "nodes": [index[id(root)]]}],
        "nodes": nodes_json,
        "meshes": meshes,
        "materials": mats_json,
        "accessors": accessors,
        "bufferViews": buffer_views,
        "buffers": [{"byteLength": len(BIN)}],
    }
    if img_view[0] is not None:
        gltf["images"] = [{"bufferView": img_view[0], "mimeType": "image/png"}]
        gltf["samplers"] = [{"magFilter": 9729, "minFilter": 9987, "wrapS": 33071, "wrapT": 33071}]
        gltf["textures"] = [{"sampler": 0, "source": 0}]
        for m in mats_json:
            if m["name"] == "NamePlate":
                m["pbrMetallicRoughness"]["baseColorTexture"] = {"index": 0}

    js = json.dumps(gltf, separators=(",", ":")).encode("utf-8")
    while len(js) % 4:
        js += b" "
    total = 12 + 8 + len(js) + 8 + len(BIN)
    with open(path, "wb") as f:
        f.write(struct.pack("<III", 0x46546C67, 2, total))
        f.write(struct.pack("<II", len(js), 0x4E4F534A))
        f.write(js)
        f.write(struct.pack("<II", len(BIN), 0x004E4942))
        f.write(BIN)
    return len(nodes_json), len(meshes)


# ---------------------------------------------------------------- preview
def collect_world(root):
    tris, mats = [], []

    def walk(n, parent):
        w = parent.dot(n.local)
        R, t = w[:3, :3], w[:3, 3]
        for (mat, geo) in n.parts:
            P, Nn, UV, Tr = geo.arrays()
            if len(P) == 0:
                continue
            tris.append((P.dot(R.T) + t)[Tr])
            mats.append(np.full(len(Tr), mat))
        for c in n.children:
            walk(c, w)
    walk(root, np.eye(4))
    return np.concatenate(tris), np.concatenate(mats)


def render(TP, TM, W=900, H=900, eye=(150, 118, 205), target=(0, 8, 0), fov=31.0,
           out="preview.png", bg=(0.90, 0.91, 0.93)):
    eye = np.array(eye, float)
    target = np.array(target, float)
    fwd = normalize(target - eye)
    right = normalize(np.cross(fwd, np.array([0, 1, 0])))
    upv = np.cross(right, fwd)
    cam = np.stack([right, upv, -fwd])

    Pc = (TP.reshape(-1, 3) - eye).dot(cam.T).reshape(-1, 3, 3)
    z = -Pc[:, :, 2]
    ok = (z > 1.0).all(axis=1)
    Pc, TM = Pc[ok], TM[ok]
    z = -Pc[:, :, 2]
    f = 1.0 / math.tan(fov * D2R / 2)
    px = ((f * Pc[:, :, 0] / z) * 0.5 + 0.5) * W
    py = (0.5 - (f * Pc[:, :, 1] / z) * 0.5) * H
    iz = 1.0 / z

    yy = np.arange(H)[:, None] / H
    img = np.array(bg)[None, None, :] * (0.94 + 0.09 * (1 - yy))[..., None]
    img = np.repeat(img, W, axis=1)
    zbuf = np.full((H, W), 1e9)

    base = np.array([[v ** (1 / 2.2) for v in m["pbr"]["baseColorFactor"][:3]] for m in MATERIALS])
    metal = np.array([m["pbr"]["metallicFactor"] for m in MATERIALS])
    rough = np.array([m["pbr"]["roughnessFactor"] for m in MATERIALS])
    L1 = normalize(np.array([-0.42, 0.78, 0.52]))
    L2 = normalize(np.array([0.80, 0.30, 0.40]))

    for ti in range(len(Pc)):
        p, vx, vy = Pc[ti], px[ti], py[ti]
        x0 = max(int(math.floor(vx.min())), 0)
        x1 = min(int(math.ceil(vx.max())), W - 1)
        y0 = max(int(math.floor(vy.min())), 0)
        y1 = min(int(math.ceil(vy.max())), H - 1)
        if x1 < x0 or y1 < y0:
            continue
        area = (vx[1] - vx[0]) * (vy[2] - vy[0]) - (vx[2] - vx[0]) * (vy[1] - vy[0])
        if abs(area) < 1e-10:
            continue
        X, Y = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
        w0 = ((vx[1] - X) * (vy[2] - Y) - (vx[2] - X) * (vy[1] - Y)) / area
        w1 = ((vx[2] - X) * (vy[0] - Y) - (vx[0] - X) * (vy[2] - Y)) / area
        w2 = 1.0 - w0 - w1
        inside = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
        if not inside.any():
            continue
        depth = w0 * iz[ti][0] + w1 * iz[ti][1] + w2 * iz[ti][2]
        sub = zbuf[y0:y1 + 1, x0:x1 + 1]
        better = inside & (depth < sub)
        if not better.any():
            continue
        sub[better] = depth[better]

        n = normalize(np.cross(p[1] - p[0], p[2] - p[0]))
        vdir = normalize(-p.mean(axis=0))
        if np.dot(n, vdir) < 0:
            n = -n
        m = TM[ti]
        cosv = abs(float(np.dot(n, vdir)))
        ndl = max(0.0, float(np.dot(n, L1))) * 0.85 + max(0.0, float(np.dot(n, L2))) * 0.35
        spec = max(0.0, float(np.dot(n, normalize(L1 + vdir)))) ** (2.0 / max(1e-3, rough[m]) ** 2)
        if metal[m] > 0.5:
            refl = 2 * float(np.dot(n, vdir)) * n - vdir
            env = 0.66 + 0.34 * max(0.0, float(refl[1])) - 0.12 * max(0.0, float(-refl[1]))
            col = base[m] * env * (0.42 + 0.80 * ndl) + np.array([1, 0.98, 0.95]) * spec * 0.9
            col += np.array([0.72, 0.78, 0.88]) * (1 - cosv) ** 4 * 0.35
        else:
            col = base[m] * (0.32 + 0.82 * ndl) + np.array([1, 1, 1]) * spec * 0.25
        img[y0:y1 + 1, x0:x1 + 1][better] = np.clip(col, 0, 1)

    out_img = Image.fromarray((np.clip(img, 0, 1) ** (1 / 2.2) * 255).astype(np.uint8))
    path = os.path.join(WORK, out)
    out_img.save(path)
    return out_img


# ---------------------------------------------------------------- main
def main():
    root = Node("ParkerManifold")
    body = Node("Body")
    build_body(body)
    label = Node("NamePlate")
    build_label(label)
    body.children.append(label)
    root.children.append(body)

    cv = Node("CenterValve", T(0, 19.0, 0))
    build_valve(cv, 1.0, M_RED, 1.0)
    root.children.append(cv)

    tilt = 34.0 * D2R
    for name, sx, ang in (("LeftAngleValve", -1, tilt), ("RightAngleValve", 1, -tilt)):
        d = np.array([math.sin(tilt) * sx, math.cos(tilt), 0.0])
        base = np.array([23.0 * sx, 13.25, 0.0]) - d * 1.6
        v = Node(name, T(base[0], base[1], base[2]).dot(Rz(ang)))
        build_valve(v, 0.90, M_BLUE, 0.95)
        root.children.append(v)

    tex, tsize = make_label_texture()
    mats = [dict(name=m["name"], pbr=dict(m["pbr"]), doubleSided=False) for m in MATERIALS]
    glb = os.path.join(OUT_DIR, "parker_manifold_3valve.glb")
    # root carries the mm -> m scale so the exported asset is metric
    n_nodes, n_meshes = write_glb(glb, root, mats, texture_png=tex,
                                  root_matrix=np.diag([0.001, 0.001, 0.001, 1.0]))

    TP, TM = collect_world(root)
    nt = len(TP)
    render(TP, TM, out="preview_front.png")
    render(TP, TM, eye=(-165, 125, 190), out="preview_left.png", fov=31.0)
    render(TP, TM, eye=(30, 235, 95), target=(0, 7, 0), out="preview_top.png", fov=31.0)
    print("GLB %s  %.1f KB  tris=%d nodes=%d meshes=%d label=%s"
          % (glb, os.path.getsize(glb) / 1024.0, nt, n_nodes, n_meshes, tsize))


if __name__ == "__main__":
    main()
