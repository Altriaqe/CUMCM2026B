"""内层布局的"可定位性"定量分析（纯几何，不用模拟器）。

对区域内每个源位置 S，问两个问题：
  1. 有几个内层点落在 S 的接收半径内（能拿到示向度）？
     接收半径取最坏值 1000m —— 用它才能保证"无论官方怎么设都成立"。
  2. 拿到的这些示向度里，最大的交会角是多少？
     交会角 >= min_angle(20°) 才可能做最小二乘定位；否则只能沿单射线取
     700m 种子点交给 homing，多付一整段横向收敛路程。

这两个数是"内层布局好不好"的直接刻画，且完全不依赖模拟器。
"""
import math
from itertools import combinations

import strategy

WORST_RECV = 1000.0
MIN_ANGLE = 20.0


def ring(r, n, phase=0.0):
    return [(r * math.cos(math.radians(phase + 360.0 * k / n)),
             r * math.sin(math.radians(phase + 360.0 * k / n)))
            for k in range(n)]


def crossing_angle(p1, a1, p2, a2):
    """两条方位线在交会点的夹角（度），取 [0,90]。"""
    d1 = math.radians(a1)
    d2 = math.radians(a2)
    # 方位差
    ang = abs(math.degrees(math.atan2(
        math.sin(d1 - d2), math.cos(d1 - d2))))
    if ang > 90.0:
        ang = 180.0 - ang
    return ang


def analyze(points, S, recv=WORST_RECV):
    """返回 (可见点数, 最大交会角)。"""
    vis = []
    for p in points:
        rho = math.hypot(S[0] - p[0], S[1] - p[1])
        if rho <= recv:
            ang = math.degrees(math.atan2(S[1] - p[1], S[0] - p[0])) % 360.0
            vis.append((p, ang))
    if len(vis) < 2:
        return len(vis), 0.0
    best = 0.0
    for (p1, a1), (p2, a2) in combinations(vis, 2):
        q = crossing_angle(p1, a1, p2, a2)
        if q > best:
            best = q
    return len(vis), best


def sample_sources(area_r=1800.0, n=3000, seed=7):
    import random
    rng = random.Random(seed)
    out = []
    while len(out) < n:
        x = rng.uniform(-area_r, area_r)
        y = rng.uniform(-area_r, area_r)
        if x * x + y * y <= area_r * area_r:
            out.append((x, y))
    return out


def report(name, points, src):
    n_vis = []
    angs = []
    unlocatable = 0
    for S in src:
        nv, q = analyze(points, S)
        n_vis.append(nv)
        angs.append(q)
        if nv >= 2 and q < MIN_ANGLE:
            unlocatable += 1
        elif nv < 2:
            unlocatable += 1
    avg_nv = sum(n_vis) / len(n_vis)
    avg_q = sum(angs) / len(angs)
    good = sum(1 for q in angs if q >= MIN_ANGLE) / len(src)
    print(f"{name:<30} 平均可见点={avg_nv:5.2f}  "
          f"平均最大交会角={avg_q:6.2f}°  "
          f"可定位(角>=20°)占比={good*100:5.1f}%  "
          f"低于阈值占比={unlocatable/len(src)*100:5.1f}%")


if __name__ == "__main__":
    src = sample_sources()
    print(f"采样 {len(src)} 个源位置，接收半径按最坏 {WORST_RECV:.0f}m，"
          f"阈值 {MIN_ANGLE}°\n")

    report("均匀 r=1123 (Kershner)", [(0.0, 0.0)] + ring(1123.0, 6), src)
    report("均匀 r=1300", [(0.0, 0.0)] + ring(1300.0, 6), src)
    report("均匀 r=1400", [(0.0, 0.0)] + ring(1400.0, 6), src)
    report("混合 400x3 + 1600x4", strategy.coverage_points_mixed()[:7], src)

    print()
    print("--- 只数'中部源'(距中心 700~1100m，接收半径内点最多的地带) ---")
    mid = [S for S in src if 700 <= math.hypot(*S) <= 1100]
    report("均匀 r=1123", [(0.0, 0.0)] + ring(1123.0, 6), mid)
    report("混合 400x3 + 1600x4", strategy.coverage_points_mixed()[:7], mid)

    print()
    print("--- 边缘源 (距中心 >1600m) ---")
    edge = [S for S in src if math.hypot(*S) > 1600]
    report("均匀 r=1123", [(0.0, 0.0)] + ring(1123.0, 6), edge)
    report("均匀 r=1300", [(0.0, 0.0)] + ring(1300.0, 6), edge)
    report("混合 400x3 + 1600x4", strategy.coverage_points_mixed()[:7], edge)
