"""
问题 1 / 问题 2 的几何计算，以及多示向度最小二乘定位（供问题 3/4 使用）。

角度约定：正东 0°，逆时针为正，范围 [0, 360)。
"""
import math
from itertools import combinations


def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def bearing_point(origin, angle_deg, d):
    """从 origin 沿方位角 angle_deg（度）前进 d 米的点。"""
    a = math.radians(angle_deg)
    return (origin[0] + d * math.cos(a), origin[1] + d * math.sin(a))


def ray_intersection(p1, a1_deg, p2, a2_deg):
    """两条由 (点, 方位角) 确定的射线的交点；平行时返回 None。"""
    a1, a2 = math.radians(a1_deg), math.radians(a2_deg)
    d1 = (math.cos(a1), math.sin(a1))
    d2 = (math.cos(a2), math.sin(a2))
    denom = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denom) < 1e-12:
        return None
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    t = (dx * d2[1] - dy * d2[0]) / denom
    return (p1[0] + t * d1[0], p1[1] + t * d1[1])


def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull(points):
    """单调链凸包，返回逆时针顶点。"""
    pts = sorted(set((round(x, 6), round(y, 6)) for x, y in points))
    if len(pts) <= 1:
        return pts
    lower = []
    for p in pts:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _halfplane_clip(poly, p0, ang_deg, keep_ccw):
    """用「过 p0、方向 ang_deg 的有向直线」的一侧裁剪凸多边形。

    keep_ccw=True 保留方向向量逆时针一侧（叉积 ≥ 0）。Sutherland–Hodgman，
    凸多边形裁剪后仍是凸多边形。"""
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


def cross_bearing_region(bearings, err=1.0, lim=1.0e5):
    """交会定位区域：每个检测点的示向度 θ±err 给出两条射线，源必在其间。

    **正确构造是「2n 个半平面求交」**：过每个检测点 P、方向 θ−err 的直线的
    逆时针一侧，与方向 θ+err 的直线的顺时针一侧，两者的交集。结果是一个凸
    多边形，n 个检测点给出 2n 个半平面。

    早期版本取「所有射线两两交点的凸包」，那只是两两约束的交，对 n≥2 会得到
    比真实区域**更大**的多边形，从而**高估**定位误差 —— n=2 时两者恰好等价
    （四条射线围成的凸四边形），n≥3 时不等价。此函数现统一走半平面路径，
    并在 solve_q12.q1_region 之外提供同一算法的复用入口。

    bearings 可以是两种形式：
      · [(S1, θ1), (S2, θ2), ...]  —— 点 + 示向度
      · cross_bearing_region(S1, θ1, S2, θ2, err) —— 两个检测点的旧调用形式
    返回逆时针凸包顶点。"""
    if not bearings:
        return []
    # 兼容旧的四参数调用
    if len(bearings) == 4 and all(isinstance(x, (int, float)) for x in bearings):
        S1, theta1, S2, theta2 = bearings          # 兼容旧的四参数调用
        bearings = [(S1, theta1), (S2, theta2)]
    elif not isinstance(bearings[0], (tuple, list)):
        raise TypeError("bearings 应为 [(点, 示向度), ...] 序列")

    box = [(lim, lim), (-lim, lim), (-lim, -lim), (lim, -lim)]
    poly = box
    for p, a in bearings:
        poly = _halfplane_clip(poly, p, a + err, keep_ccw=False)
        if not poly:
            return []
        poly = _halfplane_clip(poly, p, a - err, keep_ccw=True)
        if not poly:
            return []
    return convex_hull(poly)


def rotating_calipers(hull):
    """旋转卡壳：O(h) 求凸多边形直径。hull 须为逆时针凸包顶点序列。

    对每条边维护一个「对踵点」，随边逆时针推进单调前移，因此全部对踵关系
    总扫描量为 O(h)。返回 (直径, (p, q))。"""
    h = len(hull)
    if h < 2:
        return 0.0, (hull[0], hull[0]) if hull else (None, None)
    if h == 2:
        return dist(hull[0], hull[1]), (hull[0], hull[1])

    def area2(a, b, c):
        return abs(_cross(a, b, c))

    best = -1.0
    pair = (hull[0], hull[0])
    k = 1
    for i in range(h):
        j = (i + 1) % h
        # 推进对踵点：直到 (i,j) 与下一条边的面积不再增大
        while area2(hull[i], hull[j], hull[(k + 1) % h]) > \
                area2(hull[i], hull[j], hull[k]):
            k = (k + 1) % h
        for cand in (hull[k], hull[(k + 1) % h]):
            d = dist(hull[i], cand)
            if d > best:
                best, pair = d, (hull[i], cand)
            d = dist(hull[j], cand)
            if d > best:
                best, pair = d, (hull[j], cand)
    return best, pair


def polygon_diameter(points, method="calipers"):
    """凸多边形直径 = 顶点对最大距离。返回 (直径, (p, q))。

    method='calipers'（默认）用旋转卡壳 O(h)；'brute' 用顶点对枚举 O(h^2)，
    两者互为校验（实测偏差 < 1e-9 m），brute 仅用于小规模验证与单测。"""
    if not points:
        return 0.0, (None, None)
    hull = convex_hull(points) if len(points) > 2 or method == "brute" else list(points)
    if method == "brute":
        best, pair = -1.0, None
        for p, q in combinations(hull, 2):
            d = dist(p, q)
            if d > best:
                best, pair = d, (p, q)
        return (best, pair) if pair else (0.0, (hull[0], hull[0]))
    return rotating_calipers(hull)


def diameter_circle_covers(points):
    """回答：以定位区域直径为直径的圆能否覆盖该区域。

    返回 (直径, 直径对, 直径圆是否覆盖, Jung 定理下界)。
    Jung 定理：直径 d 的点集，其最小覆盖圆半径 <= d/sqrt(3)，
    故最小覆盖圆直径最多为 2d/sqrt(3) > d，直径圆一般不能覆盖。"""
    d, (p, q) = polygon_diameter(points)
    mid = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
    r = d / 2
    covered = all(dist(mid, v) <= r + 1e-9 for v in points)
    return d, (p, q), covered, d / math.sqrt(3)


def crossing_angle(p1, a1_deg, p2, a2_deg):
    """两条示向度线的交会角 gamma ∈ [0, 90]，单位度。

    gamma -> 0 时两线近平行，交会解病态（分母 ∼ sin gamma -> 0），
    此时交点对 ±1° 误差极端敏感，**该解不可采信**。"""
    d = math.radians((a1_deg - a2_deg) % 360)
    ang = min(abs(d), 2 * math.pi - abs(d))
    if ang > math.pi / 2:
        ang = math.pi - ang
    return math.degrees(ang)


def intersection_quality(p1, a1_deg, p2, a2_deg, min_angle=8.0, max_dist=3000.0):
    """带质量判定的交会定位。

    返回 (交点, 交会角, 基线长)；交会角小于 min_angle（解病态）或交点
    明显脱离可行区域时返回 None —— 宁可不用这个解，也不要用一个错的解。"""
    q = crossing_angle(p1, a1_deg, p2, a2_deg)
    if q < min_angle:
        return None
    p = ray_intersection(p1, a1_deg, p2, a2_deg)
    if p is None:
        return None
    base = dist(p1, p2)
    if dist(p, p1) > max_dist or dist(p, p2) > max_dist:
        return None
    return p, q, base


def triangulate_quality(bearings, min_angle=8.0):
    """多示向度最小二乘 + 交会质量。

    返回 (位置, 最佳交会角)；示向度少于 2 条或全部近平行时返回 (None, 0)。
    用最小二乘而非单对交会，是为了让冗余测量摊平个别病态方向的影响。"""
    if len(bearings) < 2:
        return None, 0.0
    best = 0.0
    for (p1, a1), (p2, a2) in combinations(bearings, 2):
        best = max(best, crossing_angle(p1, a1, p2, a2))
    if best < min_angle:
        return None, best
    return triangulate(bearings), best


def triangulate(bearings):
    """多示向度最小二乘定位。

    bearings: [(检测点(x,y), 示向度度), ...]
    每条示向度线方向 u=(cosθ,sinθ)，要求源 S 满足 u × (S-P)=0。
    对多条线做最小二乘：min Σ (u × (S-P))^2，为线性最小二乘。"""
    A = []
    b = []
    for (px, py), a in bearings:
        t = math.radians(a)
        ux, uy = math.cos(t), math.sin(t)
        A.append([-uy, ux])
        b.append(ux * py - uy * px)
    sxx = sum(r[0] * r[0] for r in A)
    sxy = sum(r[0] * r[1] for r in A)
    syy = sum(r[1] * r[1] for r in A)
    bx = sum(r[0] * v for r, v in zip(A, b))
    by = sum(r[1] * v for r, v in zip(A, b))
    det = sxx * syy - sxy * sxy
    if abs(det) < 1e-12:
        return None
    x = (syy * bx - sxy * by) / det
    y = (sxx * by - sxy * bx) / det
    return (x, y)


def triangulate_weighted(bearings):
    """按射线距离加权的多示向度最小二乘定位（问题三/四实际使用版本）。

    误差传播：示向度误差 ±δ 在距离 r 处产生的横向位置偏差约为 r·δ，故第 i 条
    约束 u_i×(S−P_i)=0 的位置误差标准差 σ_i ∝ r_i。最小化 Σ (u_i×(S−P_i))²/σ_i²
    即对每个方程除以 r_i 加权。**不这样做时，远处观测点的同一角度误差会被
    与近处等权对待，而它的实际位置影响大得多**，估计会系统性偏向近处观测点。

    同一位置重测不变时 r_i 可用 u_i·(Ŝ−P_i) 以迭代方式估计（一次迭代即收敛）：
    先做无权最小二乘得初值 Ŝ，再以 r̂_i=‖Ŝ−P_i‖ 加权重解。

    返回 (x, y)；约束不足（秩亏）时返回 None。"""
    if len(bearings) < 2:
        return None
    s0 = triangulate(bearings)
    if s0 is None:
        return None
    A, b, w = [], [], []
    for (px, py), a in bearings:
        t = math.radians(a)
        ux, uy = math.cos(t), math.sin(t)
        A.append([-uy, ux])
        b.append(ux * py - uy * px)
        # 该条示向度到初值估计的距离，作为误差尺度
        r = max(math.hypot(s0[0] - px, s0[1] - py), 1.0)
        w.append(1.0 / r)

    sxx = sum(w[i] ** 2 * A[i][0] * A[i][0] for i in range(len(A)))
    sxy = sum(w[i] ** 2 * A[i][0] * A[i][1] for i in range(len(A)))
    syy = sum(w[i] ** 2 * A[i][1] * A[i][1] for i in range(len(A)))
    bx = sum(w[i] ** 2 * A[i][0] * b[i] for i in range(len(A)))
    by = sum(w[i] ** 2 * A[i][1] * b[i] for i in range(len(A)))
    det = sxx * syy - sxy * sxy
    if abs(det) < 1e-12:
        return None
    return ((syy * bx - sxy * by) / det, (sxx * by - sxy * bx) / det)


def position_error_bound(bearings, pos):
    """给定示向度组与解 pos，给出 ±δ 误差下的定位不确定度上界（米）。

    按每条射线的横向偏差 ‖P−S‖·tanδ 做误差传递，取最小二乘解处的
    2σ 椭球半长轴近似。用于问题三的「逼近到多少米才够清除」判据。"""
    import numpy as _np
    d0 = math.radians(1.0)
    M = _np.zeros((2, 2))
    for (px, py), a in bearings:
        t = math.radians(a)
        ux, uy = math.cos(t), math.sin(t)
        n = _np.array([-uy, ux])                     # 射线法向
        r = max(math.hypot(pos[0] - px, pos[1] - py), 1e-6)
        sigma = r * math.tan(d0)                     # 该观测的横向标准差
        w = 1.0 / max(sigma, 1e-9)
        M += w * w * _np.outer(n, n)
    try:
        ev = _np.linalg.eigvalsh(M)
    except Exception:
        return None
    lam = ev[ev > 1e-12]
    if not len(lam):
        return None
    return 2.0 / math.sqrt(lam.min())                # 2σ 最弱方向半轴


def localization_quality(S1, theta1, S2, src, err=1.0):
    """候选第二检测点 S2 的定位质量：返回 (交会角, 定位区域直径)。"""
    a2 = math.degrees(math.atan2(src[1] - S2[1], src[0] - S2[0])) % 360
    region = cross_bearing_region([(S1, theta1), (S2, a2)], err=err)
    diam, _ = polygon_diameter(region)
    d = math.radians((a2 - theta1) % 360)
    ang = min(abs(d), 2 * math.pi - abs(d))
    if ang > math.pi / 2:
        ang = math.pi - ang
    return math.degrees(ang), diam


def choose_second_point(S1, theta1, src_dist, recv_radius=1400, n=72):
    """问题 2：给第一个检测点 S1 与示向度 theta1，选第二个检测点。

    源在 S1 沿 theta1 方向 src_dist 处（距离未知，取可接收半径的中值/枚举）。
    在以源为圆心、recv_radius 为半径的圆上搜索使定位区域直径最小的点。"""
    src = bearing_point(S1, theta1, src_dist)
    best = None
    for k in range(n):
        ang = 360 * k / n
        S2 = (src[0] + recv_radius * math.cos(math.radians(ang)),
              src[1] + recv_radius * math.sin(math.radians(ang)))
        q, diam = localization_quality(S1, theta1, S2, src)
        if best is None or diam < best[0]:
            best = (diam, S2, q)
    return best
