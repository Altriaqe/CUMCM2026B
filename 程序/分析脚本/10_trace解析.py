# -*- coding: utf-8 -*-
"""离线复盘：把 trace_*.jsonl 按阶段切开，定位 fail 发生在 homing 的哪一段。"""
import json, sys, collections, io, math

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

path = sys.argv[1] if len(sys.argv) > 1 else 'trace_q3_0911_224532.jsonl'
rows = [json.loads(l) for l in open(path, encoding='utf-8')]

acts = [r for r in rows if r['kind'] == 'act']

# ---- 切分：mark 只标出 homing 的「起点」，/enter 与 /exit 是天然边界 ----
# 主扫描 = 第一次 homing mark 之前的全部动作（含 /enter）；
# 每个 homing 块 = 自己的 mark 到下一个 mark 之间。
# 批处理（逐点扫）与归航之间还有一个明显信号：位置突变。
marks = [(i, r['tag']) for i, r in enumerate(rows)
         if r['kind'] == 'mark' and (r['tag'] or '').startswith('homing')]
idx_homing = [i for i, r in enumerate(rows) if r['kind'] == 'mark'
              and (r['tag'] or '').startswith('homing')]
first_hm = idx_homing[0] if idx_homing else len(rows)
sweep = [r for r in rows[:first_hm] if r['kind'] == 'act']

blocks = collections.OrderedDict()
for k, (i, tag) in enumerate(marks):
    j = marks[k + 1][0] if k + 1 < len(marks) else len(rows)
    blocks[tag] = [r for r in rows[i:j] if r['kind'] == 'act']

print('=== 主扫描 ===')
print('  动作 %d 个, vt %.1f -> %.1f' % (len(sweep), sweep[0]['vt'], sweep[-1]['vt']))
print(' ', dict(collections.Counter(r['res'] for r in sweep)))
print()

# ---- 每个频道的 homing ----
print('=== 各频道 homing ===')
print('%-40s %9s %9s %8s %6s %5s %5s %-5s %s' %
      ('频道', '起vt', '止vt', '时长', '路程', 'meas', 'clear', '结果', 'clear返回值'))
tot_homing = 0.0
for tag, rs in blocks.items():
    if not rs:
        continue
    cl = [r for r in rs if r['path'] == '/clear']
    ok = any(r['res'] == 'success' for r in cl)
    dt = rs[-1]['vt'] - rs[0]['vt']
    dist = 0.0
    for a, b in zip(rs, rs[1:]):
        if a['pos'] and b['pos']:
            dist += math.hypot(b['pos'][0] - a['pos'][0], b['pos'][1] - a['pos'][1])
    tot_homing += dt
    print('%-40s %9.1f %9.1f %8.1f %6.0f %5d %5d %-5s %s' % (
        tag[:40], rs[0]['vt'], rs[-1]['vt'], dt, dist,
        sum(1 for r in rs if r['path'] == '/measure'), len(cl),
        'OK' if ok else 'FAIL',
        dict(collections.Counter(r['res'] for r in cl))))
print()
print('homing 总时长 %.1f, 占总虚拟时间 %.1f%%' % (tot_homing, 100 * tot_homing / rows[-1]['vt']))

# ---- 核心指标：clear 动作 > 频道数 说明有频道经过了 homing 精修 ----
nch = len(set(r['ch'] for r in acts if r['ch']))
ncl = sum(1 for r in acts if r['path'] == '/clear')
print()
print('不同频道数 %d, clear 动作数 %d, homing 精修 %d 次' % (nch, ncl, ncl - nch))
print('一次性清除率 %.1f%%' % (100.0 * (2 * nch - ncl) / nch))

# ---- 展开 FAIL 的频道，逐步打印 ----
print()
for tag, rs in blocks.items():
    cl = [r for r in rs if r['path'] == '/clear']
    if any(r['res'] == 'success' for r in cl):
        continue
    print('=' * 100)
    print('FAIL:', tag)
    prev_pos = None
    for r in rs:
        step = 0.0
        if prev_pos and r['pos']:
            step = math.hypot(r['pos'][0] - prev_pos[0], r['pos'][1] - prev_pos[1])
        prev_pos = r['pos']
        print('  vt=%8.1f %-9s ch=%-3s res=%-18s pos=(%7.1f,%7.1f)'
              '  测点=(%s,%s) svd=%-6s 步长=%6.1f' % (
                  r['vt'], r['path'], r['ch'], r['res'],
                  r['pos'][0] if r['pos'] else 0, r['pos'][1] if r['pos'] else 0,
                  ('%.0f' % r['x']) if r['x'] is not None else '-',
                  ('%.0f' % r['y']) if r['y'] is not None else '-',
                  ('%.2f' % r['svd']) if r['svd'] is not None else '-',
                  step))
