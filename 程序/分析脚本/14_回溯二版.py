# -*- coding: utf-8 -*-
"""折返代价的两种口径：
   口径 A（几何）  : 每个源 P→E vs P→S，看"就地清"这条路本身值不值。
   口径 B（可比）  : 扫描结束后的实际总路程  vs  同样这 5 个点的最短路程。
"""
import json, sys, math, itertools, collections, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

rows = [json.loads(l) for l in open(sys.argv[1], encoding='utf-8')]
acts = [r for r in rows if r['kind'] == 'act']
def d(a, b): return math.hypot(a[0] - b[0], a[1] - b[1])

for i, r in enumerate(rows):
    if r['kind'] == 'mark' and (r['tag'] or '').startswith('homing'):
        first_hm = i
        break
split = sum(1 for r in rows[:first_hm] if r['kind'] == 'act')
S = tuple(acts[split - 1]['pos'])
vt_S = acts[split - 1]['vt']
vt_end = acts[-1]['vt']

# 目标点：每个频道的估计位置 E（就是 homing 实际去的点）
est = collections.OrderedDict()
for r in rows:
    if r['kind'] == 'mark' and (r['tag'] or '').startswith('homing'):
        ch = int(r['tag'].split(':ch')[1].split(':')[0])
        body = r['tag'].split('est=(')[1].split(')')[0]
        x, y = body.split(',')
        est[ch] = (float(x), float(y))

# 实际：split 之后的全部路程
actual = sum(d(acts[i]['pos'], acts[i + 1]['pos'])
             for i in range(split - 1, len(acts) - 1)
             if acts[i]['pos'] and acts[i + 1]['pos'])

# 理论：从 S 出发遍历 5 个 E 的最短路程（TSP，n=5 枚举即可）
pts = list(est.values())
best = None
for perm in itertools.permutations(range(len(pts))):
    D = d(S, pts[perm[0]]) + sum(d(pts[perm[i]], pts[perm[i + 1]])
                                 for i in range(len(perm) - 1))
    if best is None or D < best:
        best = D

print('扫描阶段结束 vt=%.1f  S=(%.0f,%.0f)' % (vt_S, S[0], S[1]))
print('整局终点       vt=%.1f' % vt_end)
print()
print('扫描后实际总路程      %8.0f m   (%.0f s)' % (actual, actual))
print('5 个估计点的最短路程  %8.0f m   (%.0f s)  ← 理论下界（含从 S 出发）' % (best, best))
print('可压缩                %8.0f m   (%.0f s，占全局 %.1f%%)'
      % (actual - best, actual - best, 100 * (actual - best) / vt_end))
print()

# 各源距离 S 的远近，决定它们是分散还是扎堆
print('各估计点到 S 的距离:')
for ch, E in est.items():
    print('  ch%-3d E=(%8.0f,%8.0f)  |S-E|=%7.0f' % (ch, E[0], E[1], d(S, E)))
print()

# 口径 A：P→E vs P→S
sight = collections.OrderedDict()
for r in acts[:split]:
    if r['path'] == '/measure' and r['res'] == 'direction' and r.get('ch') and r['pos']:
        sight.setdefault(r['ch'], []).append(tuple(r['pos']))
print('%-4s %-19s %-19s %7s %7s %7s  %s' %
      ('频道', '最后测点 P', '估计点 E', 'P→E', 'P→S', 'S→E', 'P→S 减 P→E'))
tot = 0.0
for ch in est:
    if ch not in sight:
        continue
    P = sight[ch][-1]
    E = est[ch]
    dPE, dPS, dSE = d(P, E), d(P, S), d(S, E)
    tot += dPS - dPE
    print('%-4d (%7.0f,%7.0f) (%7.0f,%7.0f) %7.0f %7.0f %7.0f  %+7.0f'
          % (ch, P[0], P[1], E[0], E[1], dPE, dPS, dSE, dPS - dPE))
print()
print('就地清逐源省下之和 %+.0f m' % tot)
