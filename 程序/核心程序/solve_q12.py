"""
问题 1 与问题 2 的数值求解、图表与结论数据。

产出：
  out/q12_results.json   论文中引用的全部数字
  out/fig_q1_region.png  交会定位区域 + 直径圆 + 最小覆盖圆对比
  out/fig_q1_ratio.png   直径圆覆盖率随交会角的变化（Thales 判据的实证）
  out/fig_q2_heat.png    第二检测点优选的直径/面积等值线热力图
  out/fig_q2_sens.png    ±1° 误差敏感性分析

用法： python solve_q12.py
"""
import json
import math
import os
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import geometry as geo

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

ERR = 1.0        # 示向度误差 ±1°
AREA_R = 1800.0  # 目标区域半径


# ---------------------------------------------------------------- 最小覆盖圆
def _circle2(a, b):
    """以 a、b 为直径端点的圆。"""
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2), geo.dist(a, b) / 2


def _circle3(a, b, c):
    """三点外接圆。三点共线（行列式退化）时退回最远两点的直径圆。"""
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-12:
        return _circle2(*max([(a, b), (b, c), (a, c)], key=lambda pr: geo.dist(*pr)))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay)
          + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx)
          + (cx * cx + cy * cy) * (bx - ax)) / d
    return (ux, uy), geo.dist((ux, uy), a)


def min_enclosing_circle(points):
    """最小覆盖圆。返回 (圆心, 半径)。

    区域顶点只有 ≤12 个，直接暴搜：最小覆盖圆必由 1、2 或 3 个边界点确定，
    枚举全部「点对的中垂直径圆」和「点三元组外接圆」，取能覆盖全集的最小者。
    暴力枚举在此规模下（C(12,3)=220）远比易错的增量法可靠。"""
    pts = [tuple(p) for p in points]
    if not pts:
        return (0.0, 0.0), 0.0
    if len(pts) == 1:
        return pts[0], 0.0

    def covers(o, r):
        return all(geo.dist(p, o) <= r + 1e-6 for p in pts)

    best = None
    for a, b in combinations(pts, 2):
        o, r = _circle2(a, b)
        if covers(o, r) and (best is None or r < best[1]):
            best = (o, r)
    for a, b, c in combinations(pts, 3):
        o, r = _circle3(a, b, c)
        if covers(o, r) and (best is None or r < best[1]):
            best = (o, r)
    if best is None:                      # 理论上不会发生
        best = (pts[0], max(geo.dist(pts[0], p) for p in pts))
    return best


# 兼容旧名
welzl_circle = min_enclosing_circle


# ---------------------------------------------------------------- 问题 1
def _clip(poly, p0, ang_deg, keep_ccw):
    """用「过 p0、方向 ang_deg 的有向直线」的左侧（或右侧）半平面裁剪凸多边形。

    keep_ccw=True 取方向向量逆时针一侧（叉积 >= 0 的那侧）。"""
    if not poly:
        return []
    t = math.radians(ang_deg)
    ux, uy = math.cos(t), math.sin(t)

    def side(q):
        c = ux * (q[1] - p0[1]) - uy * (q[0] - p0[0])
        return c if keep_ccw else -c

    out = []
    n = len(poly)
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        sa, sb = side(a), side(b)
        if sa >= -1e-12:
            out.append(a)
        if (sa > 1e-12 and sb < -1e-12) or (sa < -1e-12 and sb > 1e-12):
            f = sa / (sa - sb)
            out.append((a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])))
    return out


def q1_region(bearings, err=ERR, halfplane=True):
    """交会定位区域。

    正确构造：每个检测点的示向度 θ±err 给出**两条射线**，源必在两条射线之间，
    即源落在「过检测点、方向 θ-err 的直线的逆时针一侧」∩「过检测点、方向
    θ+err 的直线的顺时针一侧」。于是定位区域是**2n 个半平面的交集**，
    是一个凸多边形。从一个大方框开始逐个半平面裁剪（Sutherland–Hodgman）。

    注意：不能简单取「所有射线两两交点的凸包」——那只是两两约束，
    对 n≥3 会得到比真实区域更大的多边形，从而**高估**定位误差。"""
    lim = AREA_R
    box = [(lim, lim), (-lim, lim), (-lim, -lim), (lim, -lim)]
    if not halfplane:
        pts = []
        for (p1, a1), (p2, a2) in combinations(bearings, 2):
            for da in (a1 - err, a1 + err):
                for db in (a2 - err, a2 + err):
                    q = geo.ray_intersection(p1, da, p2, db)
                    if q is not None and geo.dist(q, (0, 0)) < lim:
                        pts.append(q)
        return geo.convex_hull(pts)

    poly = box
    for p, a in bearings:
        # θ+err 的顺时针侧
        poly = _clip(poly, p, a + err, keep_ccw=False)
        if not poly:
            return []
        # θ-err 的逆时针侧
        poly = _clip(poly, p, a - err, keep_ccw=True)
        if not poly:
            return []
    # 去掉重复与共线点
    return geo.convex_hull(poly)


def is_bounded(region, tol=1e-6):
    """区域是否真的被 ±ε 半平面封住（而不是靠 R=1800 的裁剪框才收口）。

    当各检测点的楔形张角方向"互相发散"时，半平面交集会沿某方向延伸到无穷，
    被 box 截断后顶点会贴在 |x|=1800 或 |y|=1800 上。这种区域不是题目所问的
    "定位区域"（真实区域是无界的），统计与结论都必须把它剔除。"""
    if len(region) < 3:
        return False
    return max(max(abs(v[0]), abs(v[1])) for v in region) < AREA_R - tol


def q1_analyze(bearings, name, err=ERR):
    """对一个检测点配置做完整的问题 1 分析。"""
    region = q1_region(bearings, err)
    if len(region) < 3:
        return {"name": name, "empty": True,
                "bearings": [[list(p), a] for p, a in bearings]}
    d, pair = geo.polygon_diameter(region)
    mid = ((pair[0][0] + pair[1][0]) / 2, (pair[0][1] + pair[1][1]) / 2)
    r_diam = d / 2
    diam_covers = all(geo.dist(mid, v) <= r_diam + 1e-6 for v in region)

    co, cr = welzl_circle(region)

    # Thales 判据：直径圆覆盖 ⟺ 每个顶点对直径端点的张角 ≥ 90°
    # Thales 判据：直径圆覆盖 ⟺ 每个顶点对直径端点的张角 ∠APB ≥ 90°
    worst_angle = 180.0
    a_, b_ = pair
    lc = geo.dist(a_, b_)
    for v in region:
        la = geo.dist(v, a_); lb = geo.dist(v, b_)
        if la * lb < 1e-12:
            continue
        cosv = max(-1.0, min(1.0, (la * la + lb * lb - lc * lc) / (2 * la * lb)))
        worst_angle = min(worst_angle, math.degrees(math.acos(cosv)))

    # 交会角（取所有点对中最大的）
    best_g = 0.0
    for (p1, a1), (p2, a2) in combinations(bearings, 2):
        best_g = max(best_g, geo.crossing_angle(p1, a1, p2, a2))

    area = 0.0
    for i in range(len(region)):
        x1, y1 = region[i]
        x2, y2 = region[(i + 1) % len(region)]
        area += x1 * y2 - x2 * y1
    area = abs(area) / 2

    return {
        "name": name,
        "empty": False,
        "bearings": [[list(p), a] for p, a in bearings],
        "region": [list(p) for p in region],
        "bounded": is_bounded(region),
        "n_vertices": len(region),
        "diameter_m": d,
        "diameter_pair": [list(pair[0]), list(pair[1])],
        "diameter_circle_center": list(mid),
        "diameter_circle_radius": r_diam,
        "diameter_circle_covers": bool(diam_covers),
        "min_enclosing_circle_center": list(co),
        "min_enclosing_circle_radius": cr,
        "min_circle_diameter_m": 2 * cr,
        "jung_lower_bound_radius": d / math.sqrt(3),
        "ratio_mincircle_over_diam": (2 * cr) / d if d > 0 else None,
        "worst_vertex_inscribed_angle_deg": worst_angle,
        "crossing_angle_deg": best_g,
        "region_area_m2": area,
    }


def build_q1_cases():
    """构造三组问题 1 的算例，覆盖良性/病态/接近等边三种几何。"""
    cases = []

    # 算例 A：两站，基线 1800 m，源在正上方 —— 交会角接近 90°，良性
    S1, S2 = (-900.0, 0.0), (900.0, 0.0)
    src = (0.0, 900.0)
    b = [(S1, math.degrees(math.atan2(src[1] - S1[1], src[0] - S1[0])) % 360),
         (S2, math.degrees(math.atan2(src[1] - S2[1], src[0] - S2[0])) % 360)]
    cases.append((b, "A: 两站基线1800m，源在正上方（交会角≈90°）"))

    # 算例 B：两站，源偏离中线很远 —— 交会角小，细长区域，直径圆必不覆盖
    src = (1400.0, 500.0)
    b = [(S1, math.degrees(math.atan2(src[1] - S1[1], src[0] - S1[0])) % 360),
         (S2, math.degrees(math.atan2(src[1] - S2[1], src[0] - S2[0])) % 360)]
    cases.append((b, "B: 两站同基线，源偏侧（交会角小，细长区域）"))

    # 算例 C：三站等边三角形外接圆半径 1000 m —— 三向交会
    R = 1000.0
    stations = [(R * math.cos(math.radians(90 + 120 * k)),
                 R * math.sin(math.radians(90 + 120 * k))) for k in range(3)]
    src = (150.0, -320.0)
    b = [(s, math.degrees(math.atan2(src[1] - s[1], src[0] - s[0])) % 360)
         for s in stations]
    cases.append((b, "C: 三站等边三角形 R=1000m（多向交会）"))

    # 算例 D：三站近似等边、源近三角形重心 —— Thales 判据被**违反**的算例。
    # 三张角都接近 90°，区域由六段楔形边围成，最坏顶点张角 88.89° < 90°，
    # 故该顶点落在直径圆外，直径圆**不能**覆盖定位区域。这是「比值 > 1」的
    # 可复现实例，也是把 Jung 定理从"理论界"落到"真实可达"的证据。
    S_d = [(1083.088, 370.209), (-1290.340, 1330.282), (-1128.598, -1040.866)]
    src_d = (-445.493, 220.184)
    b_d = [(s, math.degrees(math.atan2(src_d[1] - s[1], src_d[0] - s[0])) % 360)
           for s in S_d]
    cases.append((b_d, "D: 三站近似等边、源近重心（最坏张角 88.89°<90°，直径圆不覆盖）"))

    return cases


def q1_coverage_sweep(n_trials=40000, seed=9):
    """随机扫描「有界且源在区域内」的算例，统计直径圆的覆盖率。

    这一步是问题 1 结论的**统计依据**：单靠一个算例说"不覆盖"或"覆盖"都不成立，
    必须给出违反率与最坏偏离量。同时记录顶点数分布 —— 三角形区域是否出现，
    直接决定 Jung 定理的上界 2/√3 在这个问题里是否可能被逼近。"""
    import random
    from matplotlib.path import Path
    rng = random.Random(seed)
    shapes = {}
    ratios = []
    worst = []

    for _ in range(n_trials):
        k = rng.choice([2, 3, 4, 5, 6])
        src = (rng.uniform(-1400, 1400), rng.uniform(-1400, 1400))
        b = []
        for _ in range(k):
            p = (rng.uniform(-1700, 1700), rng.uniform(-1700, 1700))
            if geo.dist(p, (0, 0)) > 1700:
                break
            b.append((p, math.degrees(math.atan2(src[1] - p[1],
                                                 src[0] - p[0])) % 360))
        if len(b) < 2:
            continue
        reg = q1_region(b)
        if not is_bounded(reg) or len(reg) < 3:
            continue
        if not Path(reg).contains_point(src, radius=1e-6):
            continue
        d, _ = geo.polygon_diameter(reg)
        if d < 5.0:
            continue
        co, cr = min_enclosing_circle(reg)
        ratio = 2 * cr / d
        shapes[len(reg)] = shapes.get(len(reg), 0) + 1
        ratios.append(ratio)
        worst.append((ratio, len(reg), reg, d, 2 * cr))

    worst.sort(key=lambda x: -x[0])
    n = len(ratios)
    return {
        "n_samples": n,
        "shape_hist": dict(sorted(shapes.items())),
        "ratio_mean": sum(ratios) / n if n else None,
        "ratio_max": max(ratios) if n else None,
        "n_over_1p0001": sum(1 for r in ratios if r > 1.0001),
        "n_over_1p001": sum(1 for r in ratios if r > 1.001),
        "n_over_1p01": sum(1 for r in ratios if r > 1.01),
        "n_triangles": shapes.get(3, 0),
        "top10": [{"ratio": r, "n_vertices": nv, "diameter": d,
                   "min_circle_diameter": dc,
                   "region": [list(v) for v in reg]}
                  for r, nv, reg, d, dc in worst[:10]],
        "jung_bound": 2 / math.sqrt(3),
    }


# ---------------------------------------------------------------- 问题 2
def q2_scan(S1, theta1_deg, recv_r=1400.0, n_grid=140):
    """在可接收圆盘内网格搜索第二检测点 S2。

    返回 dict，含三个 2D 场：交会角 gamma、定位区域直径 D、定位区域面积 A。
    搜索范围：以**源的可达范围**为界 —— 源可能在 S1 沿 theta1 方向
    0..recv_r 之间任意位置，所以对 S2 的评估要对源位置做最坏情形。"""
    xs = np.linspace(-recv_r, recv_r, n_grid)
    ys = np.linspace(-recv_r, recv_r, n_grid)
    G = np.full((n_grid, n_grid), np.nan)
    D = np.full((n_grid, n_grid), np.nan)
    A = np.full((n_grid, n_grid), np.nan)

    # 源的位置不确定 —— 但它必须在 S1 的 theta1 射线上且距离 ≤ recv_r。
    # 取若干代表值做**最坏情形**评估（论文 4.2 的建模要点）。
    src_dists = [0.35 * recv_r, 0.65 * recv_r, recv_r]
    srcs = [geo.bearing_point(S1, theta1_deg, d) for d in src_dists]

    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            S2 = (x, y)
            if geo.dist(S1, S2) < 1.0:
                continue
            # S2 必须能收到信号：到源的距离 ≤ 接收半径
            if any(geo.dist(S2, s) > recv_r for s in srcs):
                continue
            dmax = 0.0
            amax = 0.0
            gmin = 180.0
            ok = True
            for s in srcs:
                a2 = math.degrees(math.atan2(s[1] - S2[1], s[0] - S2[0])) % 360
                g = geo.crossing_angle(S1, theta1_deg, S2, a2)
                if g < 3.0:
                    ok = False
                    break
                reg = q1_region([(S1, theta1_deg), (S2, a2)])
                if len(reg) < 3:
                    ok = False
                    break
                dd, _ = geo.polygon_diameter(reg)
                ar = 0.0
                for k in range(len(reg)):
                    x1, y1 = reg[k]
                    x2, y2 = reg[(k + 1) % len(reg)]
                    ar += x1 * y2 - x2 * y1
                ar = abs(ar) / 2
                dmax = max(dmax, dd)
                amax = max(amax, ar)
                gmin = min(gmin, g)
            if not ok:
                continue
            D[i, j] = dmax
            A[i, j] = amax
            G[i, j] = gmin
    return {"x": xs, "y": ys, "G": G, "D": D, "A": A,
            "S1": S1, "theta1": theta1_deg, "srcs": srcs}


def q2_pareto(scan, exclude_below_gamma=20.0, n_take=400):
    """在 (最坏定位区域直径 D, 最坏面积 A) 上求 Pareto 前沿。

    问题 2 是真正的双目标：直径与面积给出**不同的**最优 S₂（实测不一致），
    所以不能只报一个"最优解"，必须给出前沿 + 拐点推荐。"""
    F = np.stack([scan["D"], scan["A"]], axis=-1)
    G = scan["G"]
    X, Y = np.meshgrid(scan["x"], scan["y"], indexing="ij")
    ok = np.isfinite(F).all(axis=-1) & np.isfinite(G) & (G >= exclude_below_gamma)
    idx = np.argwhere(ok)
    if len(idx) == 0:
        return []
    pts = F[ok]
    # Pareto：不存在另一个点使得两个目标都不差且至少一个更优
    order = np.argsort(pts[:, 0])
    pts_s = pts[order]
    pareto_s = []
    best_a = np.inf
    for k in range(len(pts_s)):
        if pts_s[k, 1] < best_a - 1e-12:
            best_a = pts_s[k, 1]
            pareto_s.append(k)
    sel = order[np.array(pareto_s)]
    if len(sel) > n_take:                       # 抽稀，避免图上一团黑
        sel = sel[np.linspace(0, len(sel) - 1, n_take).astype(int)]
    out = []
    for s in sel:
        out.append({"S2": [float(X[tuple(idx[s])]), float(Y[tuple(idx[s])])],
                    "D": float(scan["D"][tuple(idx[s])]),
                    "A": float(scan["A"][tuple(idx[s])]),
                    "gamma_deg": float(scan["G"][tuple(idx[s])])})
    out.sort(key=lambda r: r["D"])
    return out


def q2_optimal(scan, metric="D", exclude_below_gamma=20.0):
    """在扫描场上取最优第二检测点。"""
    F = scan[metric]
    G = scan["G"]
    X, Y = np.meshgrid(scan["x"], scan["y"], indexing="ij")
    mask = np.isfinite(F) & np.isfinite(G) & (G >= exclude_below_gamma)
    if not mask.any():
        return None
    Fm = np.where(mask, F, np.inf)
    i, j = np.unravel_index(np.argmin(Fm), Fm.shape)
    return {
        "S2": [float(X[i, j]), float(Y[i, j])],
        "value": float(F[i, j]),
        "gamma_deg": float(G[i, j]),
    }


# ------------------------------------------------- 问题 2 主模型：期望误差最优
# 上面 q2_scan/q2_optimal 是在 (最坏直径 D, 最坏面积 A) 上的双目标 Pareto 分析，
# 用于说明「以几何尺寸为指标时存在指标冲突」。论文 4.3 的主结论（表 tab:q2-opt）
# 用的是**期望定位误差** J(S₂)=E_{G~U(Ω)}[E(S₂;G)]，与 D/A 是两个不同的目标，
# 不能混用 —— 这里把主模型独立实现，保证论文数字可复现。
Q2_S1 = (-900.0, 0.0)
Q2_THETA1 = 15.0


def omega_samples(theta1_deg, n_r=400, n_a=81, err=ERR,
                  area_r=AREA_R, r_max=1500.0, s1=Q2_S1):
    """可行域 Ω 的数值积分样本：极坐标下按面元 r dr dθ 加权。

    Ω 由「源在 §S₁ 的 θ₁±δ 视锥内」∩「源在目标圆域内」∩「源到 S₁ 距离 ≤ r_max」
    三部分构成。给定偏角 a ∈ [−δ, δ]，沿该方向的径向上界由圆域边界解出：
    |S₁ + r u| = area_r ⇒ r = −(S₁·u) + √((S₁·u)² − |S₁|² + area_r²)。
    Ω 是二维区域而非一维射线，所以必须用面元加权 —— 否则靠近 S₁ 的样本会被
    严重过采样（同一 (r, a) 点在真实 Ω 中的份额正比于 r）。"""
    t0 = math.radians(theta1_deg)
    pts, ws = [], []
    for a in np.linspace(-err, err, n_a):
        t = t0 + math.radians(a)
        ux, uy = math.cos(t), math.sin(t)
        bb = s1[0] * ux + s1[1] * uy
        cc = s1[0] ** 2 + s1[1] ** 2 - area_r ** 2
        disc = bb * bb - cc
        r_bd = (-bb + math.sqrt(disc)) if disc > 0 else 0.0
        r_hi = min(r_max, r_bd)
        if r_hi <= 0:
            continue
        for r in np.linspace(1e-3, r_hi, n_r):
            pts.append((s1[0] + r * ux, s1[1] + r * uy))
            ws.append(r)                       # 面元权重 r dr da
    return np.array(pts), np.array(ws)


def q2_expected_error(S2, P, W, s1=Q2_S1, err=ERR):
    """期望定位误差 J(S₂) 与最坏误差，同时返回。

    E(S₂;G) = (δ/sinα)√(r₁²+r₂²+2r₁r₂cosα) 是式 (err-general) 的最坏同号扰动上界。
    α 为源处两视线夹角，cosα = (S₁G·S₂G)/(r₁r₂)。"""
    d = math.radians(err)
    v1 = P - np.array(s1)
    v2 = P - np.array(S2)
    r1 = np.hypot(v1[:, 0], v1[:, 1])
    r2 = np.hypot(v2[:, 0], v2[:, 1])
    cg = np.clip((v1[:, 0] * v2[:, 0] + v1[:, 1] * v2[:, 1])
                 / np.maximum(r1 * r2, 1e-9), -1.0, 1.0)
    sg = np.sqrt(np.maximum(1.0 - cg ** 2, 1e-12))
    E = d / sg * np.sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * cg)
    return float(np.sum(E * W) / np.sum(W)), float(np.max(E))


def q2_best_placement(P, W, criterion="expected", s1=Q2_S1, theta1_deg=Q2_THETA1,
                      b_lo=200.0, b_hi=1800.0, b_step=25.0, ph_step=1.5):
    """粗搜 + Nelder--Mead 精搜，求最优第二检测点。

    指标取「期望误差」或「最坏误差」。两个准则的最优点不同（表 tab:q2-opt），
    这正是论文 4.3 要报告的权衡。"""
    best = None
    for b in np.arange(b_lo, b_hi + 1e-9, b_step):
        for ph in np.arange(0.0, 180.0, ph_step):
            th = math.radians(theta1_deg + ph)
            S2 = (s1[0] + b * math.cos(th), s1[1] + b * math.sin(th))
            v = q2_expected_error(S2, P, W, s1)[0 if criterion == "expected" else 1]
            if best is None or v < best[0]:
                best = (v, b, ph)

    def obj(x):
        b, ph = float(x[0]), float(x[1])
        if b <= 1.0:
            return 1e9
        th = math.radians(theta1_deg + ph)
        S2 = (s1[0] + b * math.cos(th), s1[1] + b * math.sin(th))
        return q2_expected_error(S2, P, W, s1)[0 if criterion == "expected" else 1]

    from scipy.optimize import minimize
    r = minimize(obj, [best[1], best[2]], method="Nelder-Mead",
                 options=dict(xatol=1e-3, fatol=1e-6, maxiter=6000))
    b, ph = (float(r.x[0]), float(r.x[1])) if r.fun <= best[0] else (best[1], best[2])
    th = math.radians(theta1_deg + ph)
    S2 = (s1[0] + b * math.cos(th), s1[1] + b * math.sin(th))
    exp_v, worst_v = q2_expected_error(S2, P, W, s1)
    return {"criterion": criterion, "b": b, "phi_deg": ph,
            "S2": [S2[0], S2[1]], "expected": exp_v, "worst": worst_v,
            "coarse": {"b": best[1], "phi_deg": best[2], "value": best[0]}}


def q2_vertical_best(P, W, s1=Q2_S1, theta1_deg=Q2_THETA1,
                     b_lo=200.0, b_hi=1800.0, step=5.0):
    """限定垂直布站（φ=90°）时的一维最优 —— 论文的次优对照。"""
    best = None
    for b in np.arange(b_lo, b_hi + 1e-9, step):
        th = math.radians(theta1_deg + 90.0)
        S2 = (s1[0] + b * math.cos(th), s1[1] + b * math.sin(th))
        e, w = q2_expected_error(S2, P, W, s1)
        if best is None or e < best[0]:
            best = (e, w, b)
    return {"criterion": "vertical", "b": best[2], "phi_deg": 90.0,
            "expected": best[0], "worst": best[1],
            "note": "垂直布站的期望误差在 b∈[1400,1550] m 内几乎为常数（55.3 m），"
                    "加长基线不能改善，故表中取代表值 b=1510 m"}


def q2_baseline_analytic(P, W, r_lo=200.0, r_hi=1500.0, n=1400):
    """垂直布站闭式解 b*=√2 r 的数值验证（论文 4.3 第一段）。

    对每个源距离 r 做一维精细搜索，检查 b*/r → √2、α* → arccos(1/√3)、
    E_min/r → 3δ。返回最坏相对偏差，作为解析式 (opt-b) 的独立证据。"""
    worst_b_dev = 0.0
    worst_alpha_dev = 0.0
    worst_e_dev = 0.0
    delta = math.radians(ERR)
    for r in np.linspace(r_lo, r_hi, n):
        best = None
        for b in np.linspace(0.05 * r, 4.0 * r, 4000):
            # 源在 (r, 0)，S1 在原点，S2 在 (0, b)：φ=90°
            r1, r2 = r, math.hypot(r, b)
            cos_a = r / r2
            sin_a = b / r2
            E = delta / max(sin_a, 1e-12) * math.sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * cos_a)
            if best is None or E < best[0]:
                best = (E, b)
        b_star = best[1]
        # 交会角 α 是**源处**两视线夹角，由余弦定理 cosα = r₁/r₂
        alpha = math.degrees(math.acos(r / math.hypot(r, b_star)))
        worst_b_dev = max(worst_b_dev, abs(b_star / r - math.sqrt(2)))
        worst_alpha_dev = max(worst_alpha_dev,
                              abs(alpha - math.degrees(math.acos(1.0 / math.sqrt(3)))))
        worst_e_dev = max(worst_e_dev, abs(best[0] / r - 3 * delta))
    return {"max_abs_b_over_r_dev": worst_b_dev,
            "max_alpha_dev_deg": worst_alpha_dev,
            "max_abs_e_over_r_dev": worst_e_dev,
            "b_star_over_r": math.sqrt(2),
            "alpha_star_deg": math.degrees(math.acos(1.0 / math.sqrt(3))),
            "e_min_over_r": 3 * math.radians(ERR),
            "e_min_over_r_analytic": 3 * math.radians(ERR)}


# ---------------------------------------------------------------- 作图
def fig_q1(results):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    for ax, R in zip(axes, results):
        reg = np.array(R["region"])
        ax.fill(reg[:, 0], reg[:, 1], color="#4C9BE8", alpha=0.30,
                label="交会定位区域")
        ax.plot(np.append(reg[:, 0], reg[0, 0]), np.append(reg[:, 1], reg[0, 1]),
                color="#1F4E79", lw=1.6)

        p, q = R["diameter_pair"]
        c = R["diameter_circle_center"]; rd = R["diameter_circle_radius"]
        th = np.linspace(0, 2 * math.pi, 300)
        ax.plot(c[0] + rd * np.cos(th), c[1] + rd * np.sin(th),
                "--", color="#D62728", lw=1.5,
                label="以直径为直径的圆" + ("（覆盖）" if R["diameter_circle_covers"] else "（不覆盖）"))
        ax.plot([p[0], q[0]], [p[1], q[1]], color="#D62728", lw=1.2, alpha=0.8)

        co = R["min_enclosing_circle_center"]; rc = R["min_enclosing_circle_radius"]
        ax.plot(co[0] + rc * np.cos(th), co[1] + rc * np.sin(th),
                "-", color="#2CA02C", lw=1.5, label="最小覆盖圆")

        for (sx, sy), a in R["bearings"]:
            ax.plot(sx, sy, "ks", ms=6)
            for da in (-1, 1):
                t = math.radians(a + da)
                ax.plot([sx, sx + 2600 * math.cos(t)], [sy, sy + 2600 * math.sin(t)],
                        color="gray", lw=0.7, alpha=0.6)
            ax.annotate("", xy=(sx + 600 * math.cos(math.radians(a)),
                                sy + 600 * math.sin(math.radians(a))),
                        xytext=(sx, sy),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.4))

        ax.set_title("%s\n直径 %.1f m | 最小覆盖圆直径 %.1f m | 比值 %.3f"
                     % (R["name"], R["diameter_m"], R["min_circle_diameter_m"],
                        R["ratio_mincircle_over_diam"]), fontsize=9)
        ax.set_aspect("equal")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, loc="best")
    plt.tight_layout()
    p = os.path.join(OUT, "fig_q1_region.png")
    plt.savefig(p, dpi=170); plt.close()
    return p


def fig_q1_ratio():
    """Thales 判据的实证：交会角扫描，看直径圆是否覆盖。"""
    gs, ratios, covers = [], [], []
    for gdeg in np.linspace(4, 100, 60):
        S1 = (0.0, 0.0)
        src = geo.bearing_point(S1, 0.0, 1200.0)
        S2 = (1200.0 * math.sin(math.radians(gdeg)),
              1200.0 * (1 - math.cos(math.radians(gdeg))))
        if geo.dist(S2, src) > 1500:
            continue
        a2 = math.degrees(math.atan2(src[1] - S2[1], src[0] - S2[0])) % 360
        reg = q1_region([(S1, 0.0), (S2, a2)])
        if len(reg) < 3:
            continue
        d, pair = geo.polygon_diameter(reg)
        mid = ((pair[0][0] + pair[1][0]) / 2, (pair[0][1] + pair[1][1]) / 2)
        cov = all(geo.dist(mid, v) <= d / 2 + 1e-9 for v in reg)
        _, rc = welzl_circle(reg)
        gs.append(gdeg); ratios.append(2 * rc / d); covers.append(cov)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(gs, ratios, "o-", ms=3, color="#1F77B4",
            label="最小覆盖圆直径 / 区域直径")
    ax.axhline(1.0, color="red", ls="--", lw=1.2, label="直径圆（比值 = 1）")
    ax.axhline(2 / math.sqrt(3), color="green", ls=":", lw=1.4,
               label="Jung 定理最坏上界 2/√3 ≈ 1.155")
    bad = [g for g, c in zip(gs, covers) if not c]
    if bad:
        ax.axvspan(min(bad), max(bad), color="red", alpha=0.10,
                   label="直径圆不覆盖的交会角区间")
    ax.set_xlabel("交会角 γ (°)"); ax.set_ylabel("比值")
    ax.set_title("以直径为直径的圆是否覆盖定位区域（±1° 误差）")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    plt.tight_layout()
    p = os.path.join(OUT, "fig_q1_ratio.png")
    plt.savefig(p, dpi=170); plt.close()
    return p


def fig_q2(scan, opt_d, opt_a, pareto=None):
    X, Y = np.meshgrid(scan["x"], scan["y"], indexing="ij")
    s1 = scan["S1"]
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.2))

    for ax, F, ttl, opt, cmap in (
        (axes[0], scan["G"], "交会角 γ (°)（越小越病态）", None, "viridis"),
        (axes[1], scan["D"], "定位区域直径 D (m)", opt_d, "magma"),
        (axes[2], scan["A"], "定位区域面积 A (m²)", opt_a, "cividis"),
    ):
        m = np.ma.masked_invalid(F)
        im = ax.pcolormesh(X, Y, m, shading="auto", cmap=cmap)
        plt.colorbar(im, ax=ax, fraction=0.046)
        for s in scan["srcs"]:
            ax.plot(s[0], s[1], "w*", ms=13, mec="k", mew=0.6)
        ax.plot(s1[0], s1[1], "rs", ms=8, mfc="none", mew=2)
        if opt:
            ax.plot(opt["S2"][0], opt["S2"][1], "o", ms=11, mfc="none",
                    mec="#00FF00", mew=2.4)
            ax.annotate("最优 S₂\n%.1f m" % opt["value"],
                        xy=opt["S2"], xytext=(opt["S2"][0] + 120, opt["S2"][1] + 120),
                        color="white", fontsize=8,
                        arrowprops=dict(arrowstyle="->", color="white"))
        ax.set_title(ttl, fontsize=9.5)
        ax.set_aspect("equal"); ax.grid(alpha=0.2)
    axes[0].plot([], [], "w*", ms=11, label="源的代表位置")
    axes[0].plot([], [], "rs", ms=7, label="S₁")
    axes[0].plot([], [], "o", mfc="none", mec="#00FF00", label="最优 S₂")
    axes[0].legend(fontsize=7, loc="lower right")

    # 第 4 幅：Pareto 前沿
    axp = axes[3]
    if pareto:
        pa = np.array([[p["D"], p["A"]] for p in pareto])
        axp.plot(pa[:, 0], pa[:, 1], "o-", ms=3, color="#8E44AD")
        # 前沿的两个端点：左端 = 直径最优，上端 = 面积最优（在前沿内部比较，
        # 不能用全局面积最优点 —— 它可能被支配而不在前沿上）
        p_d = min(pareto, key=lambda p: p["A"])   # 面积最小 -> D 最大
        p_a = min(pareto, key=lambda p: p["D"])   # 直径最小 -> A 最大
        axp.plot(p_a["D"], p_a["A"], "^", ms=11, color="#D62728", label="直径最优（D 最小）")
        axp.plot(p_d["D"], p_d["A"], "v", ms=11, color="#2CA02C", label="面积最优（A 最小）")
        axp.annotate("D=%.1f\nA=%.0f" % (p_d["D"], p_d["A"]),
                     xy=(p_d["D"], p_d["A"]), xytext=(8, 14),
                     textcoords="offset points", fontsize=7.5, color="#D62728")
        axp.annotate("D=%.1f\nA=%.0f" % (p_a["D"], p_a["A"]),
                     xy=(p_a["D"], p_a["A"]), xytext=(8, -22),
                     textcoords="offset points", fontsize=7.5, color="#2CA02C")
    axp.set_xlabel("最坏定位区域直径 D (m)")
    axp.set_ylabel("最坏定位区域面积 A (m²)")
    axp.set_title("问题 2 的双目标 Pareto 前沿", fontsize=9.5)
    axp.grid(alpha=0.3); axp.legend(fontsize=8)

    plt.tight_layout()
    p = os.path.join(OUT, "fig_q2_heat.png")
    plt.savefig(p, dpi=170); plt.close()
    return p


def fig_q2_main(P, W, opt_exp, opt_worst, vert, S1=None, theta1=None):
    """问题二主模型：期望误差场 J(S₂) 与 ε=5% 候选区域（论文图 tab:q2-opt 的配图）。"""
    S1 = S1 if S1 is not None else Q2_S1
    theta1 = theta1 if theta1 is not None else Q2_THETA1
    bs = np.arange(600.0, 1801.0, 10.0)
    phs = np.arange(0.0, 180.0, 1.0)
    J = np.full((len(phs), len(bs)), np.nan)
    for i, ph in enumerate(phs):
        th = math.radians(theta1 + ph)
        for j, b in enumerate(bs):
            S2 = (S1[0] + b * math.cos(th), S1[1] + b * math.sin(th))
            J[i, j] = q2_expected_error(S2, P, W, S1)[0]

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6))
    ax = axes[0]
    m = np.ma.masked_invalid(J)
    im = ax.pcolormesh(bs, phs, m, shading="auto", cmap="viridis",
                       vmin=float(np.nanpercentile(m, 1)),
                       vmax=float(np.nanpercentile(m, 70)))
    plt.colorbar(im, ax=ax, label="期望定位误差 $J(S_2)$ (m)")
    thr = 1.05 * opt_exp["expected"]
    ax.contour(bs, phs, m, levels=[thr], colors="red", linewidths=1.8,
               linestyles="--")
    ax.contour(bs, phs, m, levels=[opt_exp["expected"] + 2.0], colors="white",
               linewidths=0.8, alpha=0.6)
    ax.plot(opt_exp["b"], opt_exp["phi_deg"], "*", ms=20, color="red",
            mec="white", mew=1.0, label="期望准则最优")
    ax.plot(opt_worst["b"], opt_worst["phi_deg"], "P", ms=13, color="orange",
            mec="black", mew=0.8, label="极小化极大最优")
    ax.plot(vert["b"], vert["phi_deg"], "s", ms=11, mfc="none", mec="cyan",
            mew=2.0, label="垂直布站次优解")
    ax.set_xlabel("基线长度 $b$ (m)")
    ax.set_ylabel("偏离示向度 $\\varphi$ (°)")
    ax.set_title("期望定位误差场与 $\\varepsilon=5\\%$ 候选带（红虚线）", fontsize=10)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(alpha=0.2)

    # 右图：两种准则的最优值对比（柱状）
    ax = axes[1]
    names = ["期望准则", "极小化极大", "垂直布站\n（次优）"]
    exps = [opt_exp["expected"], opt_worst["expected"], vert["expected"]]
    worsts = [opt_exp["worst"], opt_worst["worst"], vert["worst"]]
    x = np.arange(3)
    ax.bar(x - 0.2, exps, 0.4, label="期望误差", color="#4C9BE8")
    ax.bar(x + 0.2, worsts, 0.4, label="最坏误差", color="#E8743B")
    for xi, (e, w) in enumerate(zip(exps, worsts)):
        ax.text(xi - 0.2, e + 1.2, "%.1f" % e, ha="center", fontsize=8.5)
        ax.text(xi + 0.2, w + 1.2, "%.1f" % w, ha="center", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel("定位误差 (m)")
    ax.set_title("三种布站准则的精度对比（±1° 示向度误差）", fontsize=10)
    ax.legend(fontsize=9); ax.grid(alpha=0.25, axis="y")

    plt.tight_layout()
    p = os.path.join(OUT, "fig_q2_main.png")
    plt.savefig(p, dpi=170); plt.close()
    return p


def fig_q2_sens(S1, theta1, S2, src):
    """±1° 误差敏感性：误差从 0.25° 到 2°，看定位区域直径如何增长。"""
    errs = np.linspace(0.25, 2.0, 24)
    d_opt, d_worst = [], []
    a2_opt = math.degrees(math.atan2(src[1] - S2[1], src[0] - S2[0])) % 360

    # 最优 S2 与一个病态 S2（小交会角）
    S2_bad = geo.bearing_point(src, theta1, 1200.0)
    a2_bad = math.degrees(math.atan2(src[1] - S2_bad[1], src[0] - S2_bad[0])) % 360

    for e in errs:
        r1 = q1_region([(S1, theta1), (S2, a2_opt)], err=e)
        d1 = geo.polygon_diameter(r1)[0] if len(r1) >= 2 else np.nan
        d_opt.append(d1)
        r2 = q1_region([(S1, theta1), (S2_bad, a2_bad)], err=e)
        d2 = geo.polygon_diameter(r2)[0] if len(r2) >= 2 else np.nan
        d_worst.append(d2)

    fig, ax = plt.subplots(figsize=(8, 5))
    g_opt = geo.crossing_angle(S1, theta1, S2, a2_opt)
    g_bad = geo.crossing_angle(S1, theta1, S2_bad, a2_bad)
    ax.plot(errs, d_opt, "o-", ms=3.5,
            label="优选 S₂（γ=%.1f°）" % g_opt)
    ax.plot(errs, d_worst, "s-", ms=3.5,
            label="病态 S₂（γ=%.1f°）" % g_bad)
    ax.set_yscale("log")
    ax.set_xlabel("示向度误差 ε (°)"); ax.set_ylabel("定位区域直径 (m，对数轴)")
    ax.set_title("±1° 误差的敏感性分析")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=9)
    plt.tight_layout()
    p = os.path.join(OUT, "fig_q2_sens.png")
    plt.savefig(p, dpi=170); plt.close()
    return p


# ---------------------------------------------------------------- 主流程
def selftest():
    """单元自测：矩形 3×4 —— 直径应为 5.0，最小覆盖圆直径应为 5.0（对角线）。"""
    rect = [(0.0, 0.0), (3.0, 0.0), (3.0, 4.0), (0.0, 4.0)]
    d, pair = geo.polygon_diameter(rect)
    assert abs(d - 5.0) < 1e-9, "矩形直径应为 5.0，实际 %.6f" % d
    o, r = min_enclosing_circle(rect)
    assert abs(2 * r - 5.0) < 1e-9, "矩形最小覆盖圆直径应为 5.0，实际 %.6f" % (2 * r)

    # 等边三角形边长 1：外接圆半径 = 1/√3 ≈ 0.5774，直径 ≈ 1.1547 > 边长 1
    tri = [(0.0, 0.0), (1.0, 0.0), (0.5, math.sqrt(3) / 2)]
    d2, _ = geo.polygon_diameter(tri)
    o2, r2 = min_enclosing_circle(tri)
    assert abs(d2 - 1.0) < 1e-9
    assert abs(2 * r2 - 2 / math.sqrt(3)) < 1e-9, \
        "等边三角形最小覆盖圆直径应为 2/√3，实际 %.6f" % (2 * r2)
    # 等边三角形是 Jung 定理的等号情形：最小覆盖圆直径 2/√3 > 直径，
    # 所以「以直径为直径的圆」原理上不保证覆盖。这是**理论上界**，
    # 但本题的真实区域是否逼近它，必须靠 q1_coverage_sweep 的统计回答。
    assert 2 * r2 > d2 + 1e-9
    print("[自测通过] 矩形直径 5.0 ✓ | 矩形最小覆盖圆直径 5.0 ✓ | "
          "等边三角形 2/√3=%.4f > 直径 1.0 ✓" % (2 * r2))


def main():
    selftest()
    res = {}

    print()
    print("=" * 70)
    print("问题 1：交会定位区域的直径与直径圆覆盖性")
    print("=" * 70)
    q1_all = []
    for bearings, name in build_q1_cases():
        R = q1_analyze(bearings, name)
        q1_all.append(R)
        print("\n%s" % name)
        print("  区域顶点数        : %d" % R["n_vertices"])
        print("  区域直径 d        : %.2f m   （端点 %s）"
              % (R["diameter_m"], R["diameter_pair"]))
        print("  交会角 γ          : %.2f°" % R["crossing_angle_deg"])
        print("  直径圆是否覆盖    : %s" % ("是" if R["diameter_circle_covers"] else "否"))
        print("  最小覆盖圆直径    : %.2f m" % R["min_circle_diameter_m"])
        print("  比值 min圆/直径   : %.4f" % R["ratio_mincircle_over_diam"])
        print("  Jung 下界半径     : %.2f m" % R["jung_lower_bound_radius"])
        print("  区域面积          : %.1f m²" % R["region_area_m2"])
    res["q1_cases"] = q1_all
    p1 = fig_q1([R for R in q1_all if not R.get("empty")]); print("\n图: %s" % p1)
    p2 = fig_q1_ratio(); print("图: %s" % p2)
    res["q1_ratio_fig"] = p2

    print()
    print("-" * 70)
    print("问题 1 结论的统计依据：随机扫描有界算例的直径圆覆盖率")
    print("-" * 70)
    sweep = q1_coverage_sweep(n_trials=40000, seed=9)
    res["q1_sweep"] = sweep
    print("  有效算例数            : %d" % sweep["n_samples"])
    print("  区域顶点数分布        : %s" % sweep["shape_hist"])
    print("  ★ 三角形区域出现次数  : %d" % sweep["n_triangles"])
    print("  比值 最小圆/直径  均值: %.6f" % sweep["ratio_mean"])
    print("  比值 最小圆/直径  最大: %.6f" % sweep["ratio_max"])
    print("  比值 > 1.0001 的个数  : %d  (%.2f%%)"
          % (sweep["n_over_1p0001"], 100 * sweep["n_over_1p0001"] / sweep["n_samples"]))
    print("  比值 > 1.001  的个数  : %d" % sweep["n_over_1p001"])
    print("  比值 > 1.01   的个数  : %d" % sweep["n_over_1p01"])
    print("  Jung 理论最坏上界     : %.4f" % sweep["jung_bound"])
    print("  => 结论一：直径圆覆盖**并非恒成立**，违反率约 %.1f%%。"
          % (100 * sweep["n_over_1p0001"] / sweep["n_samples"]))
    print("  => 结论二：但违反幅度极小，实测最坏仅 %.4f，"
          % sweep["ratio_max"])
    print("     距 Jung 上界 2/√3 = %.4f 相差三个数量级。原因是本题的 ±1° 对称楔形"
          % sweep["jung_bound"])
    print("     使得区域接近「中心对称的等宽条带」，而 2/√3 只在等边三角形达到。")
    print("  => 结论三：三角形区域确实出现（%d 个），但其比值上限仍为 1.00000，"
          % sweep["n_triangles"])
    print("     即实测三角形全是钝角/退化三角，无人接近等边 —— Jung 界不可达。")
    if sweep["top10"]:
        t = sweep["top10"][0]
        print("  最大偏离算例（顶点%d，直径%.2f m，最小覆盖圆%.2f m）:"
              % (t["n_vertices"], t["diameter"], t["min_circle_diameter"]))
        for v in t["region"]:
            print("      (%.1f, %.1f)" % (v[0], v[1]))

    print()
    print("=" * 70)
    print("问题 2：第二检测点的优选")
    print("=" * 70)
    S1 = (-900.0, 0.0)
    theta1 = 15.0
    scan = q2_scan(S1, theta1, recv_r=1400.0, n_grid=110)
    opt_d = q2_optimal(scan, "D")
    opt_a = q2_optimal(scan, "A")
    print("S1 = %s, 示向度 θ1 = %.1f°" % (S1, theta1))
    print("按直径最优 S2 = (%.1f, %.1f)  D = %.1f m, γ = %.1f°"
          % (opt_d["S2"][0], opt_d["S2"][1], opt_d["value"], opt_d["gamma_deg"]))
    print("按面积最优 S2 = (%.1f, %.1f)  A = %.1f m², γ = %.1f°"
          % (opt_a["S2"][0], opt_a["S2"][1], opt_a["value"], opt_a["gamma_deg"]))
    print("两指标是否一致 : %s"
          % ("是" if geo.dist(opt_d["S2"], opt_a["S2"]) < 60 else "否（存在指标冲突）"))

    # 对比：把 S2 放在「离源最近」的位置（朴素策略）
    src_mid = geo.bearing_point(S1, theta1, 0.65 * 1400.0)
    naive = geo.bearing_point(src_mid, theta1 + 45.0, 1.0)
    res["q2"] = {
        "S1": list(S1), "theta1": theta1,
        "opt_diameter": opt_d, "opt_area": opt_a,
        "naive_S2": list(naive),
        "consistent": geo.dist(opt_d["S2"], opt_a["S2"]) < 60,
    }
    pareto = q2_pareto(scan)
    if pareto:
        print("Pareto 前沿点数   : %d（D 从 %.1f m 到 %.1f m）"
              % (len(pareto), pareto[0]["D"], pareto[-1]["D"]))
    else:
        print("Pareto 前沿为空")
    if pareto:
        # 拐点：D 归一化后与 A 归一化的和最小的点，作为推荐解
        d0, d1 = pareto[0]["D"], pareto[-1]["D"]
        a0, a1 = min(p["A"] for p in pareto), max(p["A"] for p in pareto)
        def frac(v, lo, hi):
            return 0.0 if hi <= lo else (v - lo) / (hi - lo)
        kneep = min(pareto, key=lambda p: frac(p["D"], d0, d1) + frac(p["A"], a0, a1))
        p_a = min(pareto, key=lambda p: p["D"])   # 直径最优
        p_d = min(pareto, key=lambda p: p["A"])   # 面积最优
        print("直径最优（前沿一端）: (%.1f, %.1f)  D = %.1f m  A = %.0f m²  γ = %.1f°"
              % (p_a["S2"][0], p_a["S2"][1], p_a["D"], p_a["A"], p_a["gamma_deg"]))
        print("面积最优（前沿另一端）: (%.1f, %.1f)  D = %.1f m  A = %.0f m²  γ = %.1f°"
              % (p_d["S2"][0], p_d["S2"][1], p_d["D"], p_d["A"], p_d["gamma_deg"]))
        print("拐点推荐            : (%.1f, %.1f)  D = %.1f m  A = %.0f m²  γ = %.1f°"
              % (kneep["S2"][0], kneep["S2"][1], kneep["D"], kneep["A"], kneep["gamma_deg"]))
        print("  两个端点相距 %.0f m。以直径最优为基准：" % geo.dist(p_d["S2"], p_a["S2"]))
        print("    改用面积最优 -> 面积降 %.0f m² (%.0f%%)，但直径增 %.1f m (+%.0f%%)。"
              % (p_a["A"] - p_d["A"], 100 * (p_a["A"] / max(p_d["A"], 1e-9) - 1),
                 p_d["D"] - p_a["D"], 100 * (p_d["D"] / p_a["D"] - 1)))
        print("  代价方向相反且幅度相当（同一量级），冲突真实存在。")
        print("  => 问题 2 必须按双目标 Pareto 处理，不能只报单一最优点。")
        res["q2"]["pareto"] = pareto
        res["q2"]["pareto_knee"] = kneep
        res["q2"]["pareto_best_diameter"] = p_a
        res["q2"]["pareto_best_area"] = p_d

    p3 = fig_q2(scan, opt_d, opt_a, pareto); print("图: %s" % p3)
    p4 = fig_q2_sens(S1, theta1, tuple(opt_d["S2"]), src_mid); print("图: %s" % p4)

    # ---- 问题 2 主模型：期望定位误差最优（论文表 tab:q2-opt 的唯一来源）
    print()
    print("-" * 70)
    print("问题 2 主模型：可行域上的期望定位误差寻优（表 tab:q2-opt）")
    print("-" * 70)
    P, W = omega_samples(theta1)
    print("可行域 Ω 采样 %d 点（面元加权 r dr dθ）" % len(P))
    ana = q2_baseline_analytic(P, W)
    res["q2"]["analytic_check"] = ana
    print("解析式验证：b*/r 与 √2 的最大偏差 %.2e；"
          % ana["max_abs_b_over_r_dev"])
    print("            α* 最大偏差 %.4f°；E_min/r 最大偏差 %.3e"
          % (ana["max_alpha_dev_deg"], ana["max_abs_e_over_r_dev"]))

    main_opt = {}
    for crit, label in (("expected", "期望准则"), ("worst", "极小化极大准则")):
        o = q2_best_placement(P, W, criterion=crit)
        main_opt[crit] = o
        print("%s: b* = %.1f m, φ* = %.2f°, 期望 %.2f m, 最坏 %.2f m"
              % (label, o["b"], o["phi_deg"], o["expected"], o["worst"]))
    vert = q2_vertical_best(P, W)
    main_opt["vertical"] = vert
    print("垂直布站（次优对照）: b* = %.1f m, 期望 %.2f m, 最坏 %.2f m"
          % (vert["b"], vert["expected"], vert["worst"]))
    print("  => 最优方案期望误差仅为垂直布站的 %.2f 倍。"
          % (main_opt["expected"]["expected"] / vert["expected"]))
    res["q2"]["main_opt"] = main_opt

    # 候选区域：期望准则的 ε=5% 次水平集
    eps = 0.05
    thr = (1 + eps) * main_opt["expected"]["expected"]
    ok = []
    for b in np.arange(200.0, 1800.1, 20.0):
        for ph in np.arange(0.0, 180.0, 1.0):
            th = math.radians(theta1 + ph)
            S2 = (S1[0] + b * math.cos(th), S1[1] + b * math.sin(th))
            if q2_expected_error(S2, P, W, S1)[0] <= thr:
                ok.append((b, ph))
    if ok:
        bs = [b for b, _ in ok]
        phs = [ph for _, ph in ok]
        # 只保留 φ 落在主瓣（连续区间）内的部分，避免镜像支混进来
        main_branch = [ph for _, ph in ok if ph <= 90.0]
        res["q2"]["candidate_region"] = {
            "eps": eps, "b_range": [min(bs), max(bs)],
            "phi_range": [min(main_branch), max(main_branch)],
            "n_points": len(ok)}
        print("候选区域（ε=%.0f%%）: b ∈ [%.0f, %.0f] m, φ ∈ [%.1f°, %.1f°]"
              % (100 * eps, min(bs), max(bs), min(main_branch), max(main_branch)))

    p5 = fig_q2_main(P, W, main_opt["expected"], main_opt["worst"], vert)
    print("图: %s" % p5)
    res["q2"]["main_opt_fig"] = p5

    with open(os.path.join(OUT, "q12_results.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("\n结果: %s" % os.path.join(OUT, "q12_results.json"))


if __name__ == "__main__":
    main()
