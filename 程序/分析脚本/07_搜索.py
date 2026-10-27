"""针对性搜索：最大化 最小覆盖圆/直径 的比值。"""
import solve_q12 as s, geometry as geo, math, random
from matplotlib.path import Path

def ratio_of(bearings, src):
    r = s.q1_region(bearings)
    if len(r) < 3:
        return None
    if max(max(abs(v[0]), abs(v[1])) for v in r) >= 1799.99:
        return None
    if not Path(r).contains_point(src, radius=1e-6):
        return None
    d, _ = geo.polygon_diameter(r)
    if d < 5:
        return None
    cc = s.min_enclosing_circle(r)
    return 2*cc[1]/d, r, d, 2*cc[1]

random.seed(1234)
best = []
for trial in range(200000):
    k = random.choice([3, 4, 5, 6, 8, 10])
    # 源与检测点都尽量靠近边界，放大非对称性
    src = (random.uniform(-1700, 1700), random.uniform(-1700, 1700))
    if geo.dist(src, (0, 0)) > 1750:
        continue
    b = []
    for _ in range(k):
        p = (random.uniform(-1790, 1790), random.uniform(-1790, 1790))
        if geo.dist(p, (0, 0)) > 1790:
            break
        b.append((p, math.degrees(math.atan2(src[1]-p[1], src[0]-p[0])) % 360))
    if len(b) < 3:
        continue
    out = ratio_of(b, src)
    if out is None:
        continue
    ratio, r, d, dc = out
    best.append((ratio, k, r, d, dc, src))

best.sort(key=lambda x: -x[0])
print("搜索 %d 个有效算例，比值 top10：" % len(best))
for ratio, k, r, d, dc, src in best[:10]:
    print("  比值 %.6f  k=%2d 直径 %8.2f 最小圆Ø %8.2f 顶点%2d 源(%.0f,%.0f)"
          % (ratio, k, d, dc, len(r), src[0], src[1]))
print()
print("比值 > 1.001 的算例数：%d" % sum(1 for x in best if x[0] > 1.001))
print("比值 > 1.01  的算例数：%d" % sum(1 for x in best if x[0] > 1.01))
print("比值 > 1.05  的算例数：%d" % sum(1 for x in best if x[0] > 1.05))
r0, k0, rr, d0, dc0, src0 = best[0]
print()
print("最大比值算例的顶点：", [(round(v[0],1), round(v[1],1)) for v in rr])
print("  直径端点：", [(round(v[0],1), round(v[1],1)) for v in geo.polygon_diameter(rr)[1]])
