"""在纯几何判据下，优化覆盖点布局，把"最大角间隔"压到 180° 以下。

判据：对区域内每个源位置 S，观测点相对 S 的方位角的最大间隔 < 180°
      <=> 不存在任何朝向能让 S 处的定向源对所有观测点都背向。

这个判据完全不依赖模拟器，是纯几何的，可以直接用来选布局。
"""
import math
import random

from _gap_check import max_gap_for_source, sample_sources


def ring(r, n, phase=0.0):
    return [(r * math.cos(math.radians(phase + 360.0 * k / n)),
             r * math.sin(math.radians(phase + 360.0 * k / n)))
            for k in range(n)]


def worst_gap(points, sources):
    return max(max_gap_for_source(S, points) for S in sources)


def bad_frac(points, sources):
    return sum(1 for S in sources
               if max_gap_for_source(S, points) >= 180.0) / len(sources)


if __name__ == "__main__":
    src = sample_sources(n=1500)   # 搜索用稍小的样本，够快

    print("=== 1. 现有 25 点网的最坏源在哪 ===")
    base = [tuple(p) for p in __import__("strategy").coverage_points()]
    wg = max_gap_for_source(base[0], base)  # placeholder
    print()

    print("=== 2. 内层 7 点：环半径 × 相位扫描（无外层）===")
    best = []
    for r in range(1000, 1801, 50):
        for ph in (0.0, 30.0):
            pts = [(0.0, 0.0)] + ring(r, 6, ph)
            best.append((worst_gap(pts, src), r, ph))
    best.sort()
    for w, r, ph in best[:6]:
        print(f"  最坏间隔={w:7.2f}°  r={r}  phase={ph}")

    print()
    print("=== 3. 内层7点 + 外环(N, r) 组合扫描 ===")
    results = []
    for r_in in (1123.0, 1300.0, 1400.0, 1500.0):
        for n_out in (12, 14, 16, 18, 20, 24):
            for r_out in (1600.0, 1700.0, 1795.0):
                pts = [(0.0, 0.0)] + ring(r_in, 6) + ring(r_out, n_out, 15.0)
                results.append((worst_gap(pts, src), bad_frac(pts, src),
                                r_in, n_out, r_out, len(pts)))
    results.sort()
    print(f"  {'最坏间隔':>9} {'坏源占比':>9} {'r_in':>6} {'n_out':>6} "
          f"{'r_out':>7} {'总点数':>6}")
    for w, bf, ri, no, ro, n in results[:12]:
        print(f"  {w:9.2f} {bf*100:8.2f}% {ri:6.0f} {no:6d} {ro:7.0f} {n:6d}")

    print()
    print("=== 4. 三层布局（近 - 中 - 外）===")
    results3 = []
    for r1, n1 in ((400.0, 3), (600.0, 4)):
        for r2, n2 in ((1123.0, 6), (1300.0, 6)):
            for r3, n3 in ((1700.0, 12), (1795.0, 12), (1795.0, 16)):
                pts = ([(0.0, 0.0)] + ring(r1, n1, 60.0) + ring(r2, n2, 0.0)
                       + ring(r3, n3, 15.0))
                results3.append((worst_gap(pts, src), bad_frac(pts, src),
                                 r1, n1, r2, n2, r3, n3, len(pts)))
    results3.sort()
    for w, bf, r1, n1, r2, n2, r3, n3, n in results3[:8]:
        print(f"  最坏间隔={w:7.2f}°  坏源占比={bf*100:5.2f}%  "
              f"布局: 内{r1:.0f}×{n1} + 中{r2:.0f}×{n2} + 外{r3:.0f}×{n3} "
              f"共{n}点")

    print()
    print("=== 5. 当前实产布局复核（用大样本）===")
    src2 = sample_sources(n=4000)
    import strategy
    cur = strategy.coverage_points()
    print(f"  当前 coverage_points() 25点: 最坏间隔="
          f"{worst_gap(cur, src2):.2f}°  坏源占比={bad_frac(cur, src2)*100:.2f}%")
