"""净账计算：各种覆盖点布局的额外路程 vs 潜在收益。

只用几何和题目给定的时间常数，不依赖模拟器。
时间模型 T = L/5 + N_switch*1 + N_measure*5 + N_ok*3 + N_fail*5
"""
import math

import strategy
from _gap_check import max_gap_for_source, sample_sources


def tour_length(points, origin=(0.0, 0.0)):
    """从 origin 出发走遍所有点的路径长度（含 2-opt + Or-opt 优化）。"""
    pts = [tuple(p) for p in points]
    return _path_len(strategy.order_points_from(origin, pts), origin)


def _path_len(seq, start):
    if not seq:
        return 0.0
    L = math.hypot(seq[0][0] - start[0], seq[0][1] - start[1])
    for a, b in zip(seq, seq[1:]):
        L += math.hypot(a[0] - b[0], a[1] - b[1])
    return L


def ring(r, n, phase=0.0):
    return [(r * math.cos(math.radians(phase + 360.0 * k / n)),
             r * math.sin(math.radians(phase + 360.0 * k / n)))
            for k in range(n)]


if __name__ == "__main__":
    src = sample_sources(n=3000)
    spd = 5.0

    base = strategy.coverage_points()
    L_base = tour_length(base)
    w0 = max(max_gap_for_source(S, base) for S in src)
    b0 = sum(1 for S in src
             if max_gap_for_source(S, base) >= 180.0) / len(src)

    print(f"{'布局':<34}{'路程m':>8}{'移动s':>8}{'Δs':>8}"
          f"{'最坏间隔':>10}{'盲源%':>8}")
    print(f"{'当前 25 点(7+18)':<34}{L_base:8.0f}{L_base/spd:8.0f}"
          f"{0:8.0f}{w0:10.2f}{b0*100:8.2f}")

    # 各候选布局的额外路程（相对于当前布局的增量）
    cands = [
        ("外环 18 -> 20 点", 20),
        ("外环 18 -> 24 点", 24),
        ("外环 18 -> 30 点", 30),
        ("外环 18 -> 36 点", 36),
    ]
    for name, n in cands:
        pts = base[:7] + ring(1795.0, n, 15.0)
        L = tour_length(pts)
        w = max(max_gap_for_source(S, pts) for S in src)
        b = sum(1 for S in src
                if max_gap_for_source(S, pts) >= 180.0) / len(src)
        print(f"{name:<34}{L:8.0f}{L/spd:8.0f}{(L-L_base)/spd:8.0f}"
              f"{w:10.2f}{b*100:8.2f}")

    # 减少覆盖面：外环砍到 12 点
    pts = base[:7] + ring(1795.0, 12, 15.0)
    L = tour_length(pts)
    w = max(max_gap_for_source(S, pts) for S in src)
    b = sum(1 for S in src if max_gap_for_source(S, pts) >= 180.0) / len(src)
    print(f"{'外环 18 -> 12 点(省路程)':<34}{L:8.0f}{L/spd:8.0f}"
          f"{(L-L_base)/spd:8.0f}{w:10.2f}{b*100:8.2f}")

    print()
    print("=== 内层布局的巡回路程对比（7 点）===")
    for name, pts in (
        ("均匀 r=1123 (Kershner)", [(0.0, 0.0)] + ring(1123.0, 6)),
        ("均匀 r=1300", [(0.0, 0.0)] + ring(1300.0, 6)),
        ("均匀 r=1400", [(0.0, 0.0)] + ring(1400.0, 6)),
        ("混合 400x3 + 1600x4",
         strategy.coverage_points_mixed()[:7]),
    ):
        L = tour_length(pts)
        w = max(max_gap_for_source(S, pts) for S in src)
        b = sum(1 for S in src
                if max_gap_for_source(S, pts) >= 180.0) / len(src)
        print(f"  {name:<24} 路程={L:7.0f}m 移动={L/spd:6.0f}s "
              f"最坏间隔={w:7.2f}° 盲源={b*100:5.2f}%")
