# -*- coding: utf-8 -*-
"""按频道重组全轨迹：切出「归航段」= 该频道最后一次示向度 到 它第一次 clear。"""
import json, sys, math, collections, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

rows = [json.loads(l) for l in open(sys.argv[1], encoding='utf-8')]
acts = [r for r in rows if r['kind'] == 'act']

# 每频道最后一次 direction 的下标（该时刻就是"最后一次能定位它"）
lastsight = {}
for i, r in enumerate(acts):
    if r['path'] == '/measure' and r['res'] == 'direction' and r.get('ch'):
        lastsight[r['ch']] = i

# 扫描/归航分界：第一次出现 >400m 的位置跳变（主扫描单步最大 1123m，
# 但那只发生在环上相邻点之间；归航第一跳一定是从环点直奔估计位置）。
def jump(a, b):
    if not (a['pos'] and b['pos']):
        return 0.0
    return math.hypot(b['pos'][0] - a['pos'][0], b['pos'][1] - a['pos'][1])

split = len(acts)
for i in range(1, len(acts)):
    if jump(acts[i - 1], acts[i]) > 400:
        split = i
        break
print('主扫描 %d 个动作，vt 0 -> %.1f；此后归航' % (split, acts[split - 1]['vt']))
print()

cleared = set()
bych = collections.OrderedDict()
for i, r in enumerate(acts):
    if r['path'] != '/clear' or not r.get('ch') or r['ch'] in cleared:
        continue
    c = r['ch']
    bych[c] = acts[lastsight.get(c, i):i + 1]
    if r['res'] == 'success':
        cleared.add(c)

print('%-5s %-19s %8s %8s %7s %5s %5s  %s' %
      ('频道', '清除位置', '起vt', '止vt', '归航时长', '路程', 'meas', 'clear 动作序列'))
n1 = nmulti = 0
tot = 0.0
for c, seg in bych.items():
    clr = [r for r in seg if r['path'] == '/clear']
    meas = [r for r in seg if r['path'] == '/measure']
    dt = seg[-1]['vt'] - seg[0]['vt']
    dist = sum(math.hypot(b['pos'][0] - a['pos'][0], b['pos'][1] - a['pos'][1])
               for a, b in zip(seg, seg[1:]) if a['pos'] and b['pos'])
    if len(meas) <= 1:
        n1 += 1
    else:
        nmulti += 1
    tot += dt
    end = seg[-1]['pos']
    print('%-5d (%8.1f,%8.1f) %8.1f %8.1f %8.1f %7.0f %5d %5d  %s' % (
        c, end[0], end[1], seg[0]['vt'], seg[-1]['vt'], dt, dist, len(meas), len(clr),
        ' -> '.join(r['res'] + ('' if r['res'] == 'success' else '*') for r in clr)))

print()
print('清除 %d 个频道' % len(bych))
print('  一步到位（估计够准，直接清）: %d 个' % n1)
print('  需二次测量精修            : %d 个' % nmulti)
print('归航段总时长 %.1f s, 均值 %.1f s/频道' % (tot, tot / max(1, len(bych))))

# 全局口径
nclr = sum(1 for r in acts if r['path'] == '/clear')
nbad = sum(1 for r in acts if r['path'] == '/clear' and r['res'] != 'success')
print()
print('全局：clear 动作 %d 次，其中未命中 %d 次（%.1f%%）—— 每次代价 2s'
      % (nclr, nbad, 100.0 * nbad / nclr))
print('      （成功 3s / no_target_in_range 或 failed 5s，故空放一次只多花 2s）')
