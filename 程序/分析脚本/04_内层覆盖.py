"""纯几何：混合内层布局是否真的覆盖整个区域（按最坏接收半径 1000m）。

之前 coverage_points_mixed 的文档只声称"外层 18 点保持不变，Q3 全向源覆盖
不退化"，但从未验证**内层自己**是否覆盖。这里直接验证。
"""
import math

import strategy


def worst_nearest(points, area_r=1800.0, step=10.0):
    mx = 0.0
    arg = None
    R = int(area_r / step)
    for i in range(-R, R + 1):
        y = i * step
        for j in range(-R, R + 1):
            x = j * step
            if x * x + y * y <= area_r * area_r:
                d = min(math.hypot(x - p[0], y - p[1]) for p in points)
                if d > mx:
                    mx, arg = d, (x, y)
    return mx, arg


if __name__ == "__main__":
    inner = strategy.coverage_points_mixed()[:7]
    m, arg = worst_nearest(inner)
    print(f"混合内层 7 点（400×3 + 1600×4）")
    print(f"  区域内最大最近距离 = {m:.1f} m  出现在 {arg}")
    print(f"  {'覆盖 OK' if m <= 1000.0 else '覆盖不足 —— 存在全向源盲区'}"
          f"（判据：<= 最坏接收半径 1000m）")

    print()
    uni = [(0.0, 0.0)]
    for k in range(6):
        a = math.radians(60 * k)
        uni.append((1122.96 * math.cos(a), 1122.96 * math.sin(a)))
    m2, arg2 = worst_nearest(uni)
    print(f"均匀内层 7 点（Kershner r=1122.96）")
    print(f"  区域内最大最近距离 = {m2:.1f} m  出现在 {arg2}")
    print(f"  {'覆盖 OK' if m2 <= 1000.0 else '覆盖不足'}")

    print()
    print("混合布局里各点到圆心的距离:")
    for i, p in enumerate(inner):
        print(f"  {i}: r={math.hypot(*p):7.1f} m  ({p[0]:8.1f}, {p[1]:8.1f})")
