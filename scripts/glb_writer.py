# -*- coding: utf-8 -*-
"""
glb_writer.py - 零依赖(仅 numpy)的程序化 GLB 写入器 + 通用几何生成器。

约定：全部在毫米空间建模，Y 向上，+Z 朝观察者；写出时根节点乘 MM_TO_M
得到米制资产。材质表用 sRGB 颜色，写入前经 lin() 转到线性。

自检：  python glb_writer.py --self-test <输出目录>
      生成 demo.glb 并校验头部/材质/法线一致性。
"""
import json
import math
import os
import struct
import sys

import numpy as np

TAU = 2.0 * math.pi
D2R = math.pi / 180.0
MM_TO_M = np.diag([0.001, 0.001, 0.001, 1.0])


# ------------------------------------------------------------------ math utils
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


def Rx_(a):
    """绕 X 轴 -90 度等常用简化，保持与 Rx 一致。"""
    return Rx(a)


def S(sx, sy=None, sz=None):
    sy = sx if sy is None else sy
    sz = sx if sz is None else sz
    return np.diag([sx, sy, sz, 1.0])


def normalize(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def frame(origin, ydir, up=(0.0, 0.0, 1.0)):
    """4x4：局部 +Y 对齐 ydir，局部 +Z 尽量贴近 up。用于把回转体摆到任意方向。"""
    y = normalize(ydir)
    up = np.asarray(up, float)
    if abs(float(np.dot(y, up))) > 0.99:
        up = np.array([1.0, 0.0, 0.0])
    x = normalize(np.cross(y, up))
    z = np.cross(x, y)
    m = np.eye(4)
    m[:3, 0], m[:3, 1], m[:3, 2], m[:3, 3] = x, y, z, np.asarray(origin, float)
    return m


def inset_poly(poly, d):
    """凸多边形（XY，逆时针）向心内缩 d，返回等长点列。"""
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
    """从 origin 沿 ang 射出的射线与多边形边界最近交点距离。"""
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


# ------------------------------------------------------------------ geo buffer
class Geo:
    """三角形+逐顶点法线（可选 UV）累积器。绕序永远与法线同向。"""

    def __init__(self):
        self.P, self.N, self.UV, self.T = [], [], [], []
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

    def smooth_quad(self, pts, nrms, outward=None, uv=None):
        p = np.asarray(pts, float)
        nn = np.asarray(nrms, float)
        order = [0, 1, 2, 3]
        ngeo = np.cross(p[1] - p[0], p[2] - p[0])
        if np.linalg.norm(ngeo) < 1e-12:
            ngeo = np.cross(p[2] - p[0], p[3] - p[0])
        # 关键：按"几何法线 vs 顶点法线"定绕序，锥面/内腔才不会判反
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

    def check(self, tol=0.0):
        """校验绕序与法线是否同向。返回违规三角形数（0 才算通过）。"""
        P, N, _, Tr = self.arrays()
        bad = 0
        for (a, b, c) in Tr:
            ngeo = np.cross(P[b] - P[a], P[c] - P[a])
            if np.linalg.norm(ngeo) < 1e-10:
                continue
            if float(np.dot(normalize(ngeo), normalize(N[a] + N[b] + N[c]))) < tol:
                bad += 1
        return bad


# ------------------------------------------------------------------ generators
def g_revolve(profile, seg=64, smooth=True, cap_top=False, cap_bottom=False, flip=False):
    """按 (r, y) 轮廓车削。flip=True 翻转法线（做内腔用）。"""
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
                # 必须是环形顺序，写成 Z 序会让每个四边形自相交、一半三角朝内
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
            a0, a1 = j / seg * TAU, (j + 1) / seg * TAU
            g.flat_tri([(0, y, 0), (r * math.cos(a0), y, r * math.sin(a0)),
                        (r * math.cos(a1), y, r * math.sin(a1))], (0, 1, 0))
    if cap_bottom:
        r, y = profile[0]
        for j in range(seg):
            a0, a1 = j / seg * TAU, (j + 1) / seg * TAU
            g.flat_tri([(0, y, 0), (r * math.cos(a0), y, r * math.sin(a0)),
                        (r * math.cos(a1), y, r * math.sin(a1))], (0, -1, 0))
    if flip:
        g.N = [-np.asarray(x) for x in g.N]
        g.T = [(a, c, b) for (a, b, c) in g.T]
    return g


def g_cyl(r, y0, y1, seg=64, chamfer_top=0.0, chamfer_bot=0.0, cap=True):
    prof = []
    prof += [(r - chamfer_bot, y0), (r, y0 + chamfer_bot)] if chamfer_bot > 0 else [(r, y0)]
    prof += [(r, y1 - chamfer_top), (r - chamfer_top, y1)] if chamfer_top > 0 else [(r, y1)]
    return g_revolve(prof, seg=seg, smooth=True, cap_top=cap, cap_bottom=cap)


def g_ngon(n, af, y0, y1, phase=0.0, chamfer_top=0.0, chamfer_bot=0.0, cap=True):
    """正 n 棱柱，af = 对边距离（螺母的对边）。"""
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
            a0, a1 = angs[k] - math.pi / n, angs[k] + math.pi / n
            p = [(r0 * math.cos(a0), y0_, r0 * math.sin(a0)),
                 (r0 * math.cos(a1), y0_, r0 * math.sin(a1)),
                 (r1 * math.cos(a1), y1_, r1 * math.sin(a1)),
                 (r1 * math.cos(a0), y1_, r1 * math.sin(a0))]
            g.flat_quad(p, (math.cos(angs[k]), 0, math.sin(angs[k])))
    if cap:
        for (rr, yy, dirn) in ((rings[-1][0], rings[-1][1], 1), (rings[0][0], rings[0][1], -1)):
            for k in range(n):
                a0, a1 = angs[k] - math.pi / n, angs[k] + math.pi / n
                g.flat_tri([(0, yy, 0), (rr * math.cos(a0), yy, rr * math.sin(a0)),
                            (rr * math.cos(a1), yy, rr * math.sin(a1))], (0, dirn, 0))
    return g


def g_fluted(r, y0, y1, flutes=18, amp=0.4, seg_per=7, taper=0.55):
    """滚花/齿环。"""
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
        prof.append((r - taper * t + amp * (0.5 + 0.5 * math.cos(TAU * rings * t)) ** 2, y0 + h * t))
    return g_revolve(prof, seg=seg, smooth=False, cap_top=cap_top, cap_bottom=cap_bottom)


def g_bore(r_mouth, r_wall, depth, rings=6, amp=0.45, taper=0.5, seg=64):
    """内螺纹沉孔：局部 +Y 为钻进方向，孔底封闭。"""
    prof = [(r_mouth, 0.0), (r_wall, 1.3)]
    steps = rings * 4
    for k in range(1, steps + 1):
        t = k / steps
        prof.append((r_wall - taper * t + amp * (0.5 + 0.5 * math.cos(TAU * rings * t)) ** 2,
                     1.3 + (depth - 1.3) * t))
    g = g_revolve(prof, seg=seg, smooth=False, flip=True)
    r_end = prof[-1][0]
    for j in range(seg):
        a0, a1 = j / seg * TAU, (j + 1) / seg * TAU
        g.flat_tri([(0, depth, 0), (r_end * math.cos(a0), depth, r_end * math.sin(a0)),
                    (r_end * math.cos(a1), depth, r_end * math.sin(a1))], (0, -1, 0))
    return g


def emit_face_hole(g, origin, u, v, quad2d, hole_c, hole_r, outward, seg=72):
    """在 (origin,u,v) 平面上的矩形里挖一个圆孔，输出环形三角带。"""
    o, u, v = np.asarray(origin, float), np.asarray(u, float), np.asarray(v, float)
    out, hc = np.asarray(outward, float), np.asarray(hole_c, float)

    def p3(pt):
        return o + u * pt[0] + v * pt[1]

    for j in range(seg):
        a0, a1 = j / seg * TAU, (j + 1) / seg * TAU
        i0 = hc + hole_r * np.array([math.cos(a0), math.sin(a0)])
        i1 = hc + hole_r * np.array([math.cos(a1), math.sin(a1)])
        o0 = hc + ray_poly_hit(quad2d, hc, a0) * np.array([math.cos(a0), math.sin(a0)])
        o1 = hc + ray_poly_hit(quad2d, hc, a1) * np.array([math.cos(a1), math.sin(a1)])
        g.flat_quad([p3(i0), p3(i1), p3(o1), p3(o0)], out)


def strip_plate(L=33.0, W=16.0, thick=6.4, holes=((-11.3, 3.05), (0.0, 2.9), (11.3, 3.05)),
                lobes=6, amp=0.085, nx=76, bev=0.85):
    """长圆板 + 若干圆孔（手轮/压板）。holes=[(x, r)]，孔边缘顶点精确落在圆上。

    沿长轴切条带，每条带扣除孔造成的上下缺口；孔壁用同样的弦端点接缝，不会裂缝。
    """
    R = W / 2.0
    st = L / 2.0 - R

    def zout(x):
        ax = abs(x)
        z = R if ax <= st else math.sqrt(max(0.0, R * R - (min(ax, L / 2.0) - st) ** 2))
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
            I0, I1 = O0 - out * bev, O1 - out * bev
            g.flat_quad([O0 + [0, ytop, 0], O1 + [0, ytop, 0],
                         I1 + [0, ytop - bev, 0], I0 + [0, ytop - bev, 0]], out)
            g.flat_quad([I0 + [0, ytop - bev, 0], I1 + [0, ytop - bev, 0],
                         I1 + [0, ybot + bev, 0], I0 + [0, ybot + bev, 0]], out)
            g.flat_quad([I0 + [0, ybot + bev, 0], I1 + [0, ybot + bev, 0],
                         O1 + [0, ybot, 0], O0 + [0, ybot, 0]], out)
    return g


def chamfered_prism(poly, zf, zb, bevel, front_hole=None, side_ports=(), seg=72, port_depth=14.0):
    """凸多边形棱柱 + 前后倒角，可选：正面一个圆孔、若干侧面圆孔(自动带沉孔)。

    返回 [(role, Geo)]，role 取 'face' | 'side' | 'cavity'，由调用方映射到材质。
    front_hole = ((x, y), r) 或 None
    side_ports = [(edge_index, t, zd, r)]：t 沿该边从首顶点算起，
                 zd 为该边深度方向(从后倒角边往前)的孔心位置。
    """
    out = []
    ins = inset_poly(poly, bevel)
    g_front = Geo()
    if front_hole is not None:
        emit_face_hole(g_front, (0, 0, zf), (1, 0, 0), (0, 1, 0), ins, front_hole[0], front_hole[1],
                       (0, 0, 1), seg=seg)
    else:
        c = np.mean(ins, axis=0)
        for i in range(len(ins)):
            a, b = ins[i], ins[(i + 1) % len(ins)]
            g_front.flat_tri([(c[0], c[1], zf), (a[0], a[1], zf), (b[0], b[1], zf)], (0, 0, 1))
    out.append(("face", g_front))

    g_back = Geo()
    c = np.mean(ins, axis=0)
    for i in range(len(ins)):
        a, b = ins[i], ins[(i + 1) % len(ins)]
        g_back.flat_tri([(c[0], c[1], zb), (a[0], a[1], zb), (b[0], b[1], zb)], (0, 0, -1))
    out.append(("face", g_back))

    ports = {p[0]: p for p in side_ports}
    for i in range(len(poly)):
        a, b = np.array(poly[i]), np.array(poly[(i + 1) % len(poly)])
        e = b - a
        L = float(np.linalg.norm(e))
        u = e / L
        outn = np.array([u[1], -u[0], 0.0])
        gb_ = Geo()
        gb_.flat_quad([(a[0], a[1], zf - bevel), (b[0], b[1], zf - bevel),
                       (ins[(i + 1) % len(ins)][0], ins[(i + 1) % len(ins)][1], zf),
                       (ins[i][0], ins[i][1], zf)], outn + np.array([0, 0, 1.0]))
        gb_.flat_quad([(b[0], b[1], zb + bevel), (a[0], a[1], zb + bevel),
                       (ins[i][0], ins[i][1], zb),
                       (ins[(i + 1) % len(ins)][0], ins[(i + 1) % len(ins)][1], zb)],
                      outn + np.array([0, 0, -1.0]))
        out.append(("face", gb_))

        gw = Geo()
        if i in ports:
            _, t, zd, r = ports[i]
            depth = (zf - bevel) - (zb + bevel)
            emit_face_hole(gw, (a[0], a[1], zb + bevel), (u[0], u[1], 0), (0, 0, 1),
                           [(0, 0), (L, 0), (L, depth), (0, depth)], (t, zd), r, outn, seg=seg)
            out.append(("side", gw))
            mouth = np.array([a[0] + u[0] * t, a[1] + u[1] * t, zb + bevel + zd])
            bore = g_bore(r, r * 0.78, port_depth, rings=6, amp=r * 0.06, taper=r * 0.06, seg=64)
            out.append(("cavity", Geo().extend(bore, frame(mouth, -outn))))
        else:
            gw.flat_quad([(a[0], a[1], zb + bevel), (b[0], b[1], zb + bevel),
                          (b[0], b[1], zf - bevel), (a[0], a[1], zf - bevel)], outn)
            out.append(("side", gw))
    return out


# ------------------------------------------------------------------ materials
def lin(c):
    """sRGB 0-1 -> glTF 需要的线性值。"""
    return [round(max(0.0, min(1.0, v)) ** 2.2, 5) for v in c]


def mat(name, srgb, metallic=0.0, roughness=0.5, tex=False):
    pbr = {"baseColorFactor": lin(srgb) + [1.0], "metallicFactor": metallic,
           "roughnessFactor": roughness}
    return {"name": name, "pbr": pbr, "tex": tex}


STEEL = mat("StainlessSteel", (0.90, 0.905, 0.915), 0.95, 0.30)
MACHINED = mat("StainlessMachined", (0.89, 0.895, 0.905), 0.95, 0.36)
CAVITY = mat("ThreadCavity", (0.30, 0.31, 0.33), 0.70, 0.50)
PAPER = mat("NamePlate", (1.0, 1.0, 1.0), 0.0, 0.55, tex=True)


# ------------------------------------------------------------------ scene tree
class Node:
    def __init__(self, name, matrix=None):
        self.name = name
        self.local = np.eye(4) if matrix is None else matrix
        self.children = []
        self.parts = []          # [(material_index, Geo)]

    def add(self, mat_index, geo, matrix=None):
        if matrix is not None:
            geo = Geo().extend(geo, matrix)
        if len(geo.P):
            self.parts.append((mat_index, geo))
        return self

    def world(self, parent=np.eye(4)):
        return parent.dot(self.local)


# ------------------------------------------------------------------ GLB writer
def write_glb(path, root, materials, texture_png=None, root_matrix=MM_TO_M):
    """写出 GLB 2.0。materials 用上面的 mat() 结构；texture_png 给铭牌贴图路径。"""
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

    nodes_json, index = [], {}
    for n in order:
        prims = []
        for (mi, geo) in n.parts:
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
            accessors.append({"bufferView": nv, "componentType": 5126, "count": int(len(Nn)),
                              "type": "VEC3"})
            attrs["NORMAL"] = len(accessors) - 1
            if geo.has_uv:
                tv = add_view(UV.astype(np.float32).tobytes(), 34962)
                accessors.append({"bufferView": tv, "componentType": 5126, "count": int(len(UV)),
                                  "type": "VEC2"})
                attrs["TEXCOORD_0"] = len(accessors) - 1
            iv = add_view(Tr.astype(np.uint32).tobytes(), 34963)
            accessors.append({"bufferView": iv, "componentType": 5125, "count": int(Tr.size),
                              "type": "SCALAR"})
            prims.append({"attributes": attrs, "indices": len(accessors) - 1,
                          "material": mi, "mode": 4})
        if prims:
            meshes.append({"name": n.name, "primitives": prims})
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

    # 键名必须是 pbrMetallicRoughness，否则查看器静默退化成白色默认材质
    mats_json = [{"name": m["name"], "pbrMetallicRoughness": dict(m["pbr"]),
                  "doubleSided": False} for m in materials]

    gltf = {
        "asset": {"version": "2.0", "generator": "photo-to-glb glb_writer"},
        "scene": 0,
        "scenes": [{"name": root.name, "nodes": [index[id(root)]]}],
        "nodes": nodes_json,
        "meshes": meshes,
        "materials": mats_json,
        "accessors": accessors,
        "bufferViews": buffer_views,
        "buffers": [{"byteLength": len(BIN)}],
    }
    if img_view[0] is not None:
        gltf["images"] = [{"bufferView": img_view[0], "mimeType": "image/png"}]
        gltf["samplers"] = [{"magFilter": 9729, "minFilter": 9987,
                             "wrapS": 33071, "wrapT": 33071}]
        gltf["textures"] = [{"sampler": 0, "source": 0}]
        for m in mats_json:
            if any(x["name"] == m["name"] and x.get("tex") for x in materials):
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


def verify_glb(path):
    """重新读回并校验：magic/版本/长度/JSON 可解析/材质键名。"""
    d = open(path, "rb").read()
    mag, ver, tot = struct.unpack("<III", d[:12])
    jl, jt = struct.unpack("<II", d[12:20])
    js = json.loads(d[20:20 + jl])
    problems = []
    if mag != 0x46546C67:
        problems.append("bad magic")
    if ver != 2:
        problems.append("bad version")
    if tot != len(d):
        problems.append("length mismatch")
    for m in js.get("materials", []):
        if "pbrMetallicRoughness" not in m:
            problems.append("material %s missing pbrMetallicRoughness" % m.get("name"))
    for me in js.get("meshes", []):
        for p in me["primitives"]:
            ai = p["attributes"].get("POSITION")
            if ai is not None and "min" not in js["accessors"][ai]:
                problems.append("POSITION accessor without min/max")
    tris = sum(js["accessors"][p["indices"]]["count"] // 3
               for me in js["meshes"] for p in me["primitives"])
    return {"bytes": len(d), "nodes": len(js["nodes"]), "meshes": len(js["meshes"]),
            "materials": [m["name"] for m in js.get("materials", [])],
            "tris": tris, "problems": problems}


# ------------------------------------------------------------------ self test
def _self_test(outdir):
    os.makedirs(outdir, exist_ok=True)
    root = Node("Demo")
    body = Node("Body")
    poly = [(-33.0, -20.0), (33.0, -20.0), (33.0, 6.5), (13.0, 20.0), (-13.0, 20.0),
            (-33.0, 6.5)]
    parts = chamfered_prism(poly, 20.0, -20.0, 1.6, front_hole=((-1.5, -8.0), 5.8),
                            side_ports=[(1, 13.0, 18.4, 7.0)])
    role_map = {"face": 1, "side": 1, "cavity": 2}
    for item in parts:
        body.add(role_map[item[0]], item[1])
    valve = Node("Valve", T(0, 19.0, 0))
    valve.add(0, g_ngon(6, 20.0, 0.0, 5.0, chamfer_top=1.1))
    valve.add(0, g_cyl(7.5, 4.0, 11.6, seg=56, chamfer_top=1.0, chamfer_bot=0.6))
    valve.add(0, g_threaded(5.3, 15.2, 18.4, rings=4, amp=0.22, seg=48))
    valve.add(3, g_fluted(6.9, 21.9, 25.4, flutes=16, amp=0.35, seg_per=6, taper=0.5))
    valve.add(0, g_cyl(2.15, 28.8, 34.2, seg=24))
    valve.add(0, strip_plate(), T(0, 34.4, 0))
    root.children += [body, valve]

    bad = sum(g.check() for n in root.children for (_, g) in n.parts)
    path = os.path.join(outdir, "demo.glb")
    nn, nm = write_glb(path, root, [STEEL, MACHINED, CAVITY, mat("Red", (0.86, 0.16, 0.11), 0.0, 0.42)])
    info = verify_glb(path)
    print("winding violations:", bad)
    print("glb:", path, info)
    assert bad == 0, "绕序/法线不一致"
    assert not info["problems"], info["problems"]
    assert info["tris"] > 1000
    print("self-test OK")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        i = sys.argv.index("--self-test")
        _self_test(sys.argv[i + 1] if len(sys.argv) > i + 1 else ".")
    else:
        print(__doc__)
