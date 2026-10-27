import solve_q12 as s, geometry as geo, math, random
from matplotlib.path import Path

BOUND = 1799.99

def analyze(seed, trials):
    random.seed(seed)
    stat = {"bounded": 0, "noncyc": 0, "bad": [], "worst": []}
    for _ in range(trials):
        k = random.choice([2, 3, 4, 5, 6, 7])
        src = (random.uniform(-1400, 1400), random.uniform(-1400, 1400))
        b = []
        for _ in range(k):
            p = (random.uniform(-1700, 1700), random.uniform(-1700, 1700))
            if geo.dist(p, (0, 0)) > 1700:
                break
            b.append((p, math.degrees(math.atan2(src[1]-p[1], src[0]-p[0])) % 360))
        if len(b) < 2:
            continue
        r = s.q1_region(b)
        if len(r) < 3:
            continue
        # 剔除被 1800 裁剪框截断
        if max(max(abs(v[0]), abs(v[1])) for v in r) >= BOUND:
            continue
        # 无绕向依赖的点在凸多边形内判定
        pth = Path(r)
        if not pth.contains_point(src, radius=1e-6):
            continue
        d, _ = geo.polygon_diameter(r)
        if d < 5:
            continue
        cc = s.min_enclosing_circle(r)
        dev = max(abs(geo.dist(v, cc[0]) - cc[1]) for v in r)
        ratio = 2*cc[1]/d
        stat["bounded"] += 1
        if dev > 1e-6:
            stat["noncyc"] += 1
            if len(stat["bad"]) < 4:
                stat["bad"].append((dev, len(r), r, d, 2*cc[1], ratio))
        if ratio > 1.0001:
            stat["worst"].append((ratio, len(r), r, d, 2*cc[1]))
    return stat

for seed in (9, 21, 33):
    st = analyze(seed, 60000)
    n = st["bounded"]
    print("seed %2d: 有界有效算例 %6d | 顶点不共圆 %6d (%.1f%%) | 直径圆不覆盖 %d (%.2f%%)"
          % (seed, n, st["noncyc"], 100*st["noncyc"]/max(n,1),
             len(st["worst"]), 100*len(st["worst"])/max(n,1)))
    st["worst"].sort(key=lambda x: -x[0])
    for ratio, nv, r, d, dc in st["worst"][:2]:
        print("    不覆盖 比值 %.5f 直径 %.2f 最小圆Ø %.2f 顶点%d" % (ratio, d, dc, nv))
        print("       顶点 %s" % [(round(v[0],1), round(v[1],1)) for v in r])
