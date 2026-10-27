"""纯几何覆盖判据检查（不依赖模拟器）。

判据（来自参考材料）：定向源检测域是 180° 半圆盘，所以对区域内任意一点 S，
若所有观测点相对 S 的方位角中，相邻两个的最大间隔 >= 180°，就存在一个朝向
使 S 处**没有任何观测点落进锥内** -> 该源永远测不到。

即：最大角间隔 < 180° 是"零漏测"的必要条件。

用法：python _gap_check.py
"""
import math
import random

import strategy


def max_gap_for_source(S, points):
    """源位于 S 时，各观测点相对 S 的方位角的最大间隔（度）。"""
    angs = sorted(math.degrees(math.atan2(p[1] - S[1], p[0] - S[0])) % 360.0
                  for p in points)
    if len(angs) < 2:
        return 360.0
    gaps = [b - a for a, b in zip(angs, angs[1:])]
    gaps.append(360.0 - angs[-1] + angs[0])
    return max(gaps)


def sample_sources(area_r=1800.0, n=4000, seed=12345):
    """区域内均匀随机采样（含靠近边界的情形）。"""
    rng = random.Random(seed)
    out = []
    while len(out) < n:
        x = rng.uniform(-area_r, area_r)
        y = rng.uniform(-area_r, area_r)
        if x * x + y * y <= area_r * area_r:
            out.append((x, y))
    return out


def evaluate(name, points, sources, area_r=1800.0):
    worst = 0.0
    worst_pt = None
    n_bad = 0
    for S in sources:
        g = max_gap_for_source(S, points)
        if g > worst:
            worst, worst_pt = g, S
        if g >= 180.0:
            n_bad += 1
    rho = math.hypot(*worst_pt) if worst_pt else 0.0
    print(f"{name:<28} 点数={len(points):>3}  "
          f"最大角间隔={worst:7.2f}°  "
          f">=180°的源占比={n_bad/len(sources)*100:5.2f}%  "
          f"最坏源距中心={rho:7.1f}m")
    return worst, worst_pt


if __name__ == "__main__":
    src = sample_sources()
    print(f"采样 {len(src)} 个源点，区域半径 1800m\n")

    inner_mixed = strategy.coverage_points_mixed()[:7]
    inner_uniform = strategy.coverage_points()[:7]
    full_mixed = strategy.coverage_points_mixed()
    full_uniform = strategy.coverage_points()

    evaluate("Q3 混合内层(7点)", inner_mixed, src)
    evaluate("Q3 均匀内层(7点)", inner_uniform, src)
    evaluate("Q3 混合+外18环(25点)", full_mixed, src)
    evaluate("Q3 均匀+外18环(25点)", full_uniform, src)

    print()
    # 参考材料的方案：7 点内层(r=1122.96) + 18 点外层(r=1795)，与我们的对比
    ref = strategy.coverage_points()[:7] + full_uniform[7:]
    evaluate("参考材料 7+18(25点)", ref, src)

    # 内层只用混合的 7 点，外层保持
    print()
    print("内层半径扫描（7 点均匀环，无外层）：")
    for r in (1000, 1122.96, 1300, 1400, 1500, 1600, 1700):
        pts = [(0.0, 0.0)] + [
            (r * math.cos(math.radians(60 * k)),
             r * math.sin(math.radians(60 * k))) for k in range(6)]
        evaluate(f"  均匀环 r={r:.0f}", pts, src)
