# -*- coding: utf-8 -*-
"""算折返代价：early_clear 门槛放行后，边缘源能否省下归航路程。

对内层只能拿到 1 条示向度的源：
  A. 现状 —— 扫描继续走完剩下的环点，从终点 (S) 折返到估计位置 E
  B. 放行 —— 在该源最后一个测点 (P) 就地转到 E 清掉，再回到 S 继续

B - A 的净收益 = dist(S,E) - dist(P,E) - dist(E,S)
             = -dist(P,E)   ← 因为从 E 回 S 两个方案一样，只差"要不要绕 P→E"
等等，这里要小心：B 方案在 P 点分叉去 E，清完必须回 S。
A 方案是 P→S 直接回，再从 S→E。
所以 B 总路程 = dist(P,E) + dist(E,S)
   A 总路程 = dist(P,S) + dist(S,E) = dist(P,S) + dist(S,E)
两者之差（B-A）= dist(P,E) + dist(E,S) - dist(P,S) - dist(S,E)
             = dist(P,E) - dist(P,S)
即：只要 dist(P,E) < dist(P,S)，就地清就划得来。
本局 S 是主扫描终点，不是内层环点，按实际轨迹取。
"""
import json, sys, math, collections, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

rows = [json.loads(l) for l in open(sys.argv[1], encoding='utf-8')]
acts = [r for r in rows if r['kind'] == 'act']


def d(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


# ---- 分界：第一个 homing 标记。此前的所有动作（内层扫 + early_clear +
#              外层补盲 + 扫描中就地的 clear）都算"扫描阶段"。
for i, r in enumerate(rows):
    if r['kind'] == 'mark' and (r['tag'] or '').startswith('homing'):
        first_hm = i
        break
else:
    first_hm = len(rows)
split = sum(1 for r in rows[:first_hm] if r['kind'] == 'act')
S = acts[split - 1]['pos']
vt_S = acts[split - 1]['vt']
print('扫描阶段结束于 vt=%.1f，位置 S = (%.1f, %.1f)   动作数=%d' % (vt_S, S[0], S[1], split))
print()

# ---- 每频道在主扫描期间的全部示向度测点 ----
sight = collections.OrderedDict()
for r in acts[:split]:
    if r['path'] == '/measure' and r['res'] == 'direction' and r.get('ch') and r['pos']:
        sight.setdefault(r['ch'], []).append((r['pos'][0], r['pos'][1], r['vt']))

cleared_in_sweep = set()
for r in acts[:split]:
    if r['path'] == '/clear' and r['res'] == 'success' and r.get('ch'):
        cleared_in_sweep.add(r['ch'])

# ---- 归航段（split 之后）的估计点 E 与实测终点 ----
homing = collections.OrderedDict()
for r in rows:
    if r['kind'] == 'mark' and (r['tag'] or '').startswith('homing'):
        homing[r['tag']] = []

est = {}
for tag in homing:
    ch = int(tag.split(':ch')[1].split(':')[0])
    body = tag.split('est=(')[1].split(')')[0]
    x, y = body.split(',')
    est[ch] = (float(x), float(y))

print('%-4s %-8s %-20s %-20s %8s %8s %8s %8s  %s' % (
    '频道', '示向度', '最后测点 P', '估计点 E', 'P→E', 'P→S', 'S→E', '实际走', '就地清净省'))
print('-' * 118)

tot_save = 0.0
for ch in sorted(est):
    pts = sight.get(ch, [])
    if not pts:
        continue
    P = (pts[-1][0], pts[-1][1])
    E = est[ch]
    dPE = d(P, E)
    dPS = d(P, S)
    dSE = d(S, E)
    idx = [i for i, r in enumerate(acts)
           if r['path'] == '/clear' and r.get('ch') == ch and r['res'] == 'success'][0]
    actual = sum(d(acts[i]['pos'], acts[i + 1]['pos'])
                 for i in range(split - 1, idx) if acts[i]['pos'] and acts[i + 1]['pos'])
    save = dPS - dPE          # 就地清省下的净路程
    tot_save += save
    print('%-4d %-8d (%7.1f,%7.1f) (%7.1f,%7.1f) %8.0f %8.0f %8.0f %8.0f  %+8.0f' % (
        ch, len(pts), P[0], P[1], E[0], E[1], dPE, dPS, dSE, actual, save))

print()
print('5 个源就地清可省净路程合计 %+.0f m  (%.0f s 虚拟时间)' % (tot_save, tot_save))
print()
print('说明：P = 该源最后一个示向度测点；S = 主扫描终点；E = homing 用的估计位置。')
print('      「就地清净省」= P→S − P→E，为正说明从 P 直奔 E 比先回 S 再折返更近。')
