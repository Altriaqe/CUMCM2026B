"""
问题 3 / 问题 4 的搜索、定位与清除策略。

核心思路（两阶段）：
1. 覆盖扫描：在能覆盖整个目标区域（按最坏接收半径 1000m 设计）的一组检测点上，
   依次检测 20 个频道，发现所有干扰源并采集多条示向度。
2. 定位 + 清除：对每个频道用最小二乘定位得到估计位置，再按近邻顺序逐个
   “测向—逼近—清除”（homing）。
"""
import math
from itertools import combinations

import geometry as geo


AREA_R = 1800.0       # 目标区域半径
WORST_RECV = 1000.0   # 最坏有效接收半径（覆盖设计用）


def coverage_points(area_r=AREA_R, worst_recv=WORST_RECV,
                    outer_r=None, outer_n=18, outer_phase=0.0,
                    inner_phase=0.0, inner_r=None):
    """返回覆盖整个圆域的检测点集合。两层布局。

    === 内层（必须）===
    圆心 + 半径 inner_r 上的 6 个点（正六边形 + 圆心 = Kershner 构型）。
    把最坏内部覆盖距离对 inner_r 做 20m 网格扫描（area_r=1800）：

        inner_r=1000.00  最坏 1059.4 m  ✗ 边缘漏 59m
        inner_r=1050.00  最坏 1033.9 m  ✗
        inner_r=1100.00  最坏 1010.2 m  ✗
        inner_r=1122.96  最坏 1000.0 m  ✓ 临界，恰好等于接收半径下界
        inner_r=1300.00  最坏  936.5 m  ✓ 覆盖更松，但环长 7800m > 6738m

    最坏点恒在六个边界极值处（0,±area_r 及其旋转），其中 (0,-area_r) 到
    最近环点 (-r/2, -r*√3/2) 的距离满足

        d² = r² - √3·area_r·r + area_r²

    令 d = worst_recv 解出 r = (√3·area_r - √(4·worst_recv² - area_r²)) / 2。
    area_r=1800、worst_recv=1000 时 r = 1122.9558 —— 这就是 1122.96 的来历，
    它是「覆盖刚好够到接收半径下界」的那个根，**不是** Kershner 半径
    （Kershner 下界是 2/√3·1000 = 1154.7，那是「一个盘盖多大区域」的问题，
    这里问的是「7 个盘盖住固定区域」，两者不同）。

    所以 1122.96 覆盖「刚好够」：再小就有边缘空隙，再大则环变长而白白多走路
    （1200 时环周长 7200m、1300 时 7800m，而 1123 时只有 6738m）。
    默认值直接取 inner_r = area_r - worst_recv/2 = 1300 是只看覆盖、不看路程
    的结果，比最优点多走 1062m —— 调用方应显式传 1122.96。
    （旧 docstring 写「最坏距离 936.5 m」对应的是 inner_r=1300，已改正。）

    === 外层（针对定向源，默认构造但 robot.py 只取前 7 点）===
    注：本函数**总是**返回 18 个外环点，所以调用方若只要内层必须切片
    `[:7]`（robot.py 正是这么做的，Q3/Q4 都不加外环 —— 实测 13 点比 7 点
    多走 3647m 而那一局外层一个频道都没需要补测）。外层的作用在于：

    定向源只覆盖 180° 锥，判据是「检测点相对源的方向与源朝向的点积 >= 0」，
    等价于以源为界的一个半平面。于是定向源的**盲区就是它背面的半圆**。
    当源靠近区域边缘、朝向又指向区域外时，整个盲区半圆落在 1800 m 之外 ——
    此时内层 7 个点**全部在盲区里**，无论怎么调整内层参数都测不到它。
    实测：13 源 / 定向比 0.4 的 10 个种子里，7 个漏掉的源全部是这种几何
    （距中心 1404~1794 m，7 个内层点命中锥内的个数 = 0）。

    外层 18 点取 r = area_r - 5，把检测点铺到**贴着边界**的地方，从「源的
    外侧」补一条视轴。机器人本身就被限制在这个圆盘内，走到 1795 不额外
    付代价，却把最恶劣的「锥精确径向外指」情形（rho>1600 且朝向±15°内）
    的可见率从 34.7% 拉到 97.8%（外环点数与半径的联合扫描，见注释末）：
        r=1750,n=12 -> 34.7%   r=1795,n=12 -> 93.6%   r=1795,n=18 -> 97.8%
    外环从 12 点加到 18 点只多走 375 m（巡回 12016 -> 12393 m），换来
    4.2 个百分点的最恶劣情形覆盖，划算。边缘源可见率（随机朝向）：62.2% -> 97.0%。

    时间代价：19 点而非 7 点，但 detection_sweep 对已集满 min_sightings
    的频道会提前跳过，实测总时长仍在预算内。"""
    if inner_r is None:
        # 闭式最优半径（见 docstring 推导）：既让最坏覆盖恰达 worst_recv，
        # 又让环最短。旧值 area_r - worst_recv/2 = 1300 比它多走 1062m。
        inner_r = (math.sqrt(3.0) * area_r
                   - math.sqrt(max(4.0 * worst_recv ** 2 - area_r ** 2, 0.0))) / 2.0
    pts = [(0.0, 0.0)]
    for k in range(6):
        a = math.radians(inner_phase + 60 * k)
        pts.append((inner_r * math.cos(a), inner_r * math.sin(a)))
    if outer_r is None:
        outer_r = area_r - 5.0
    for k in range(outer_n):
        a = math.radians(outer_phase + 360.0 * k / outer_n)
        pts.append((outer_r * math.cos(a), outer_r * math.sin(a)))
    return pts


def coverage_points_mixed(area_r=AREA_R, outer_n=18, outer_phase=0.0,
                          near_r=1122.96, mid_r=1122.96, phase=0.0):
    """Q3 覆盖点集合：内层 7 点 + 外层若干点。

    === 为什么内层回到「圆心 + 6 点均匀环 r=1122.96」===
    曾经的「近圈 3×r=400 + 中远圈 4×r=1600」混合布局经纯几何复核后被否掉：

    1. **有覆盖洞。** 混合内层到区域内的最大最近距离是 1312 m（出现在
       (1270, -1270) 方向），超过最坏接收半径 1000 m —— 那一带有全向源
       会漏。均匀 r=1122.96 的最大最近距离恰为 1000.0 m，是 7 点覆盖的
       最小可行半径。
    2. **边缘源反而更差。** 源在 r=1700/方位20°、接收半径 1250 m 时，
       混合内层只剩 1 个可见点（交会角 0°，无法定位），均匀 r=1122.96
       有 2 个可见点、交会角 71.5°。混合布局把点摊到两圈后，「第二个
       观测点够不够得着」这件事先崩了。
    3. **更贵。** 内层巡回路程 9357 m vs 6738 m，白多走 2619 m（+523 s）。

    均匀环半径：最坏覆盖距离随 r **单调递减**（r=900→1115m，r=1559→900m
    为极小点），所以「覆盖到 1000m」不是给出一个区间，而是给出一个**下界**
    r >= 1122.96。环半径越大覆盖越松、但环也越长，于是最省路程的选择就是
    刚刚满足覆盖的那个 r = 1122.96（闭式推导见 coverage_points docstring，
    它不是 2/√3×1000 —— 那个 1154.7 是另一个问题的最优值）。

    === 外层 ===
    对**全向源**（Q3）内层已独立覆盖全区，外层只用于兜住「内层拿不到
    第二条示向度」的频道；对**定向源**（Q4）外层是补盲主力，必须够密。
    调用方按题目传 outer_n。"""
    if near_r != mid_r:                     # 兼容旧调用：混合布局已废弃
        raise ValueError("混合内层布局已废弃（存在覆盖洞），请用均匀环")
    inner_r = near_r
    pts = [(0.0, 0.0)]
    for k in range(6):
        a = math.radians(phase + 60.0 * k)
        pts.append((inner_r * math.cos(a), inner_r * math.sin(a)))
    outer_r = area_r - 5.0
    for k in range(outer_n):
        a = math.radians(outer_phase + 360.0 * k / outer_n)
        pts.append((outer_r * math.cos(a), outer_r * math.sin(a)))
    return pts


def coverage_points_grid(step=600, extent=2400):
    """Q4 专用覆盖点：(2*extent/step+1)² 个均匀方格点。

    数学保证：区域内任意一点到最近网格点不超过 step*√2/2 ≈ 424m < 1000m
    （最坏接收半径下界）。对任意朝向的定向源，其 600m 网格单元的四个顶点
    凸包包含至少一个落在源的正面半平面内的点，从而保证至少一次 direction
    反馈。step=600, extent=2400 → 81 个点，路程约 52000m（蛇形排布后）。"""
    pts = []
    xs = list(range(-extent, extent + 1, step))
    for i, y in enumerate(range(-extent, extent + 1, step)):
        row = xs if i % 2 == 0 else xs[::-1]   # 蛇形减少移动
        for x in row:
            pts.append((float(x), float(y)))
    return pts


def verify_coverage(points, area_r=AREA_R, worst_recv=WORST_RECV, step=50):
    """采样验证：区域圆内每个采样点是否都在某检测点 worst_recv 范围内。"""
    bad = 0
    total = 0
    y = -area_r
    while y <= area_r:
        x = -area_r
        while x <= area_r:
            if x * x + y * y <= area_r * area_r:
                total += 1
                if min(geo.dist((x, y), p) for p in points) > worst_recv:
                    bad += 1
            x += step
        y += step
    return total, bad


def detection_sweep(client, points, channels=range(1, 21), min_sightings=2,
                    auto_clear=False, cleared=None, no_signal_limit=None,
                    no_signal_min_r=0.0, divert=False,
                    divert_min_r=950.0, divert_min_save=150.0,
                    divert_seed_dist=700.0):
    """在每个检测点按频道顺序检测，收集每个频道的 (检测点, 示向度) 列表。

    频道一旦集满 min_sightings 条有效示向度就跳过（省切换+检测时间）。

    auto_clear: True 时碰到 near 立刻 clear，省去后续 homing 整段。
    cleared:    传入一个 set，auto_clear 成功的频道写入其中；None 则忽略。
    no_signal_limit: 单频道累计 no_signal 达到此次数后跳过该频道（用于 Q3
        全向源；Q4 定向源请传 None，因为定向源朝外时合法地给出大量 no_signal）。
    no_signal_min_r: 仅当检测点距原点 >= 此半径时才应用 no_signal_limit。
        Q3 内层近圈只有 r=400，边缘源（r≈1700）距近圈点约 1300m，超出
        1000m 接收半径，那里的 no_signal 是合法的、不代表空频道。中远圈
        r=1600 的点距边缘源仅 200-500m，连续无信号才真正指向空频道。

    divert: True 时启用「顺路分流」—— 实验性，**实测净负，robot.py 默认关**。
        动机是好的：内层环 r=1122.96 到边缘源（r≈1700）只有 700m，源却在环
        外，主扫描走完后机器人往往停在**对侧**，要穿场 1200~3200m 回去清它
        （实测一局 5 个边缘源合计穿场 10469m）。但实现方式错了，见下。

    === 为什么顺路分流净亏（同种子 A/B，14 源）===
      种子 42 -6.1%  种子 21 -5.4%  种子 13 -0.3%  种子 7 +1.3%  种子 99 +24.1%
    均值净负，且尾部风险极大（seed 99 多走 6622m = +30%）。

    两个原因叠加：
    (1) 判据算错对象。divert_worth 把「环上必须走的那一段」当成省下的路程，
        于是**源越远越倾向分流**，方向反了。真正的纯绕行差只有 271/337m 量级
        （种子 99 的 ch20/ch6），而检测点位置本身有 ±150m 噪声 —— 不可分辨。
    (2) 分流会**打断内层环的连续遍历**。7 点是按 Kershner 半径排的环，
        相邻 1123m；一旦在环点 i 处出环，回到环上时落在 i+1 或更远，
        环遍历本身从 6738m 涨到 14413m。种子 99 的 OFF 轨迹是干净的
        1123×6 环，ON 轨迹变成 1686/2337/2274/2937/2570 的折返段。

    结论：**要省的是「扫描结束后从 S 出发的穿场」，不在检测点处出环**
    （出环的代价 ≥ 穿场）。若将来重做，正确方向是把「扫描路径 + 清除顺序」
    当成一个联合 TSP 求解（清除点作为可选节点插进环），而不是逐点贪心分流。

    触发条件（三条全部满足才分流）：
      1. 该频道当前只有 1 条示向度 —— 多条时走 early_clear 的主路径；
      2. 该示向度所在检测点距原点 >= divert_min_r（默认 950m）；
      3. divert_worth(...) > divert_min_save。

    返回 sightings: dict[channel] = [(point, svd_deg) or (point, None), ...]
    None 表示 near 命中（距离 <=5m）。"""
    if cleared is None:
        cleared = set()
    sightings = {c: [] for c in channels}
    no_sig = {c: 0 for c in channels}
    npts = len(points)
    for i, (x, y) in enumerate(points):
        apply_limit = (no_signal_limit is not None and
                       math.hypot(x, y) >= no_signal_min_r)
        for c in channels:
            if c in cleared:
                continue
            if apply_limit and no_sig[c] >= no_signal_limit:
                continue
            if sum(1 for s in sightings[c] if s[1] is not None) >= min_sightings:
                continue
            r = client.measure(x, y, c)
            mr = r.get("measure_result")
            if mr == "direction":
                sightings[c].append(((x, y), r["svd_deg"]))
                no_sig[c] = 0
            elif mr == "near":
                sightings[c].append(((x, y), None))
                if auto_clear:
                    if client.clear(x, y, c).get("clear_result") == "success":
                        cleared.add(c)
            else:
                no_sig[c] += 1

        # === 本检测点扫完，做一次顺路分流机会评估 ===
        if not divert or i + 1 >= npts:
            continue
        p_now = getattr(client, "position", (x, y)) or (x, y)
        p_wait = points[i + 1]
        cands = []
        for c in channels:
            if c in cleared or not sightings[c]:
                continue
            valid = [s for s in sightings[c] if s[1] is not None]
            # 条件 1：只有 1 条示向度 —— 多条时留给 early_clear 的
            # 交会角门槛处理，那里比单条射线可靠得多。
            if len(valid) != 1:
                continue
            # 条件 2：该示向度取自边缘带检测点
            if math.hypot(valid[0][0][0], valid[0][0][1]) < divert_min_r:
                continue
            seed = single_bearing_seed(valid, seed_dist=divert_seed_dist)
            if seed is None:
                continue
            save = divert_worth(p_now, seed, p_wait)
            if save > divert_min_save:
                cands.append((save, c, seed))
        if not cands:
            continue
        # 多个频道同时够格时只做省得最多的那一个：清一个就要走一趟，
        # 其余频道在后续检测点会重新评估（位置更新，判据也更准）。
        cands.sort(reverse=True)
        _, ch, seed = cands[0]
        init_h = [(p, a) for (p, a) in sightings[ch] if a is not None]
        if homing_clear(client, ch, seed, init_hist=init_h,
                        max_probes=12):
            cleared.add(ch)
    return sightings


def order_points_from(origin, points, use_two_opt=True):
    """从 origin 出发、贪心最近邻 + 2-opt + Or-opt(1/2/3) 优化的开放路径重排。

    移动是本题时间的第一大头（5 m/s），路径质量直接决定虚拟时间。"""
    if not points:
        return []

    order = []
    remaining = list(points)
    cur = origin
    while remaining:
        i = min(range(len(remaining)), key=lambda k: geo.dist(cur, remaining[k]))
        order.append(remaining.pop(i))
        cur = order[-1]

    if not use_two_opt or len(order) < 3:
        return order

    def plen(seq, start):
        if not seq:
            return 0.0
        L = geo.dist(start, seq[0])
        for a, b in zip(seq, seq[1:]):
            L += geo.dist(a, b)
        return L

    best, best_len = list(order), plen(order, origin)

    def _2opt(seq, seq_len):
        improved = True
        while improved:
            improved = False
            for i in range(len(seq) - 1):
                for j in range(i + 1, len(seq)):
                    cand = seq[:i] + seq[i:j + 1][::-1] + seq[j + 1:]
                    L = plen(cand, origin)
                    if L < seq_len - 1e-9:
                        seq, seq_len, improved = cand, L, True
        return seq, seq_len

    def _oropt(seq, seq_len, seg_len):
        improved = True
        while improved:
            improved = False
            n = len(seq)
            for i in range(n - seg_len + 1):
                seg = seq[i:i + seg_len]
                rest = seq[:i] + seq[i + seg_len:]
                for j in range(len(rest) + 1):
                    if j in (i, i + 1):
                        continue
                    cand = rest[:j] + seg + rest[j:]
                    L = plen(cand, origin)
                    if L < seq_len - 1e-9:
                        seq, seq_len, improved = cand, L, True
                        break
                if improved:
                    break
        return seq, seq_len

    for _ in range(10):
        prev_len = best_len
        best, best_len = _2opt(best, best_len)
        best, best_len = _oropt(best, best_len, 1)
        best, best_len = _oropt(best, best_len, 2)
        best, best_len = _oropt(best, best_len, 3)
        if best_len >= prev_len - 1e-9:
            break

    return best


def interleaved_sweep_clear(client, points, channels=range(1, 21),
                            min_sightings=2, inner_n=7, max_delay=3,
                            min_angle=25.0):
    """边扫边清：按扫描点序走，一旦某频道的定位置信度够高就**就地清除**。

    === 为什么这是最大的一块收益 ===
    时间模型 T = 移动/5 + 切频道 + 5×检测 + 5×清除。实测扫描段 5300 s 里
    约 4000 s 是移动，而「先扫完 25 个点、再统一折返清除」意味着每个源都要
    从扫描路径外**额外往返一次**。边扫边清把清除的起点直接设在当前扫描点上：
    源通常就在附近（扫到它的时候机器人离它 ≤1000 m），省掉整段折返。

    === 触发条件（不能一有两条示向度就清）===
    两条示向度在交会角过小时解是病态的，过早触发只会让 homing 白跑。
    要求交会角 >= min_angle；此外对**只有 1 条示向度**的频道，只有当它已经
    连续 max_delay 个扫描点都只拿到这一条（说明后续点测不到它，多半是定向源
    锥外）时才提前处理 —— 否则继续扫，多攒一条示向度比急着清更划算。

    === 已清除频道永久剔除 ===
    清除成功的频道不再参与任何后续测量（题设：每频道最多 1 个源）。

    返回 (sightings, cleared)：sightings 是全部测量记录（供复盘），
    cleared 是成功清除的频道集合。"""
    sightings = {c: [] for c in channels}
    cleared = set()
    alive = [c for c in channels]
    last_seen_at = {c: -99 for c in channels}   # 该频道最后一次拿到示向度的点序
    pending = {}                                # 待清除的频道 -> 估计位置
    origin = getattr(client, "position", (0.0, 0.0)) or (0.0, 0.0)
    order = order_points_from(origin, points)

    for idx, (x, y) in enumerate(order):
        for c in list(alive):
            r = client.measure(x, y, c)
            mr = r.get("measure_result")
            if mr == "near":
                # 已在 5m 内且锥内：就地清除，必成
                if client.clear(x, y, c).get("clear_result") == "success":
                    cleared.add(c)
                    alive.remove(c)
                    pending.pop(c, None)
                continue
            if mr != "direction":
                continue
            sightings[c].append(((x, y), r["svd_deg"]))
            last_seen_at[c] = idx

        # 判断哪些频道该立即处理
        for c in list(alive):
            ss = sightings[c]
            valid = [(p, a) for (p, a) in ss if a is not None]
            if not valid:
                continue
            enough = len(valid) >= min_sightings
            # 单条示向度且连续 max_delay 个点没有新增 -> 提前处理（多半是定向盲区）
            stale_single = (len(valid) < min_sightings
                            and idx - last_seen_at[c] >= max_delay)
            if not (enough or stale_single):
                continue
            pos = geo.triangulate(valid) if len(valid) >= 2 else None
            if pos is None:
                if len(valid) == 1:
                    pos = geo.bearing_point(valid[0][0], valid[0][1], 700.0)
                else:
                    continue
            if len(valid) >= 2:
                q = 0.0
                for (p1, a1), (p2, a2) in combinations(valid, 2):
                    q = max(q, geo.crossing_angle(p1, a1, p2, a2))
                if q < min_angle:
                    continue
            # 就地清除：从当前位置出发
            if homing_clear(client, c, pos, min_angle=min_angle):
                cleared.add(c)
                alive.remove(c)
                pending.pop(c, None)
                last_seen_at.pop(c, None)
        # 位置归位：homing 会把机器人带走，扫描点序以点名而非当前位置驱动，
        # 但 client.position 已被 homing 改写，下一步的移动距离由模拟器按
        # 真实位置算，无需在此手动纠正。

    return sightings, cleared


def batched_sweep_clear(client, points, channels=range(1, 21),
                        min_sightings=2, inner_n=7, min_angle=25.0):
    """分批清扫：内环扫完先清一轮，再扫外环、补清。

    === 与「边扫边清」的区别（实测教训）===
    逐点触发 homing 反而比先扫后清**慢 10%**（638.9 vs 578.6 s/源），因为
    每个频道在「第一次凑够示向度」的那一刻就被处理，homing 的起点由扫描
    进度而非源的位置决定，常常停在区域的另一侧；而先扫后清时清除顺序是
    2-opt 全局最优的。串行化把 TSP 的收益吃掉了。

    分批则两者兼得：内环只有 7 个点、路程短，扫完后一次性用 2-opt 排出
    清除顺序（此时能覆盖 88% 的源），清完再去外环补盲。清除开头的折返只
    发生两次而不是每个源一次。

    返回 (sightings, cleared)。"""
    sightings = {c: [] for c in channels}
    cleared = set()
    inner, outer = split_scan_points(points, inner_n)
    origin = getattr(client, "position", (0.0, 0.0)) or (0.0, 0.0)

    for group in (inner, outer):
        if not group:
            continue
        start = getattr(client, "position", origin) or origin
        alive = [c for c in channels if c not in cleared]
        if not alive:
            break
        group = order_points_from(start, group)
        for (x, y) in group:
            for c in list(alive):
                r = client.measure(x, y, c)
                mr = r.get("measure_result")
                if mr == "near":
                    if client.clear(x, y, c).get("clear_result") == "success":
                        cleared.add(c)
                        alive.remove(c)
                    continue
                if mr == "direction":
                    sightings[c].append(((x, y), r["svd_deg"]))
        # 本批扫完：对已发现的频道做一次全局排序清除
        est = {}
        for c in alive:
            valid = [(p, a) for (p, a) in sightings[c] if a is not None]
            if not valid:
                continue
            pos = geo.triangulate(valid) if len(valid) >= 2 else None
            if pos is None or (len(valid) >= 2 and
                               max(geo.crossing_angle(p1, a1, p2, a2)
                                   for (p1, a1), (p2, a2) in combinations(valid, 2))
                               < min_angle):
                if len(valid) == 1:
                    pos = geo.bearing_point(valid[0][0], valid[0][1], 700.0)
                elif pos is None:
                    continue
            if pos is not None:
                est[c] = pos
        here = getattr(client, "position", origin) or origin
        for ch, x, y in clear_order(here, est):
            init_h = [(p, a) for (p, a) in sightings[ch] if a is not None]
            if homing_clear(client, ch, (x, y), min_angle=min_angle,
                            init_hist=init_h):
                cleared.add(ch)
                if ch in alive:
                    alive.remove(ch)
    return sightings, cleared


def split_scan_points(points, inner_n=7):
    """把覆盖点分成「内层 / 外层」两组，外层只用来补盲。"""
    return points[:inner_n], points[inner_n:]


def divert_worth(p_now, p_src, p_wait, min_save=0.0):
    """[判据已被实测否掉，保留供参考；detection_sweep(divert=True) 默认不再启用]

    原意：判断「此刻就在 p_now 处分流去清 p_src」是否比「先走到 p_wait 再回来」近。

    它写成：

        省下 = |p_now-p_wait| + |p_wait-p_src| - |p_now-p_src|

    **这个式子算错了对象。** |p_now-p_wait| 是环上必须走的一段（机器人无论如何
    都要依次路过每个检测点），不是分流能省下的东西。真正的代价是**纯绕行差**：

        绕行差 = |p_now-p_src| + |p_src-p_wait| - |p_now-p_wait|

    两者判断方向常常相反：上式在源离环越远时越大，于是越远的源越容易被判成
    「值得分流」，而越远的源恰恰是绕行差最大的那个。实测种子 99 一局，
    ch20 绕行差 271m、ch6 337m —— 而检测点位置本身有 ±150m 噪声，这个量级
    根本不可分辨，判据却把它们全判成「省 >150m」。

    返回上式的结果（米）。"""
    return (geo.dist(p_now, p_wait) + geo.dist(p_wait, p_src)
            - geo.dist(p_now, p_src))


def single_bearing_seed(hist, seed_dist=700.0, area_r=AREA_R):
    """只有 1 条示向度时的源位置种子点：沿该方位角射线走 seed_dist。

    取 hist 里**最近一条**示向度（最新测量，误差最小）。示向度带 ±1° 误差，
    在 700m 处横向偏差约 12m —— 远小于 20m 清除半径之外的余量，所以这个
    种子点足以让 homing 一两次测量内收敛（实测 16 个源里 5 个这样处理，
    全部一次成功，没有一次走 recover 散布序列）。"""
    if not hist:
        return None
    (p, a) = hist[-1]
    if a is None:
        return None
    return geo.bearing_point(p, a, seed_dist)


def _range_along_bearing(p, ang_deg, area_r):
    """估计源距：从 p 沿 ang 的射线与区域边界圆的交点距离。

    这是源距的**上界**（源在区域内），作为先验足够——补测点只需要量级正确。
    射线完全背离边界时退回 area_r / 2。"""
    # 圆 |p + t*u| = area_r 的最小正根
    ux = math.cos(math.radians(ang_deg))
    uy = math.sin(math.radians(ang_deg))
    b = p[0] * ux + p[1] * uy
    c = p[0] * p[0] + p[1] * p[1] - area_r * area_r
    disc = b * b - c
    if disc < 0.0:
        return area_r / 2.0
    t = -b + math.sqrt(disc)
    return t if t > 1.0 else area_r / 2.0


def _bearing_to(p, q):
    """p 看向 q 的方位角（度）。"""
    return math.degrees(math.atan2(q[1] - p[1], q[0] - p[0])) % 360.0


def _bearing_estimate(valid):
    """用已有的 >=2 条示向度估计源位置。

    取所有示向度对的交点的中位数（对异常值稳健），比最小二乘简单且够用。
    返回 None 表示无法估计（示向度对几乎全平行）。"""
    xs, ys = [], []
    for (p1, a1), (p2, a2) in combinations(valid, 2):
        if geo.crossing_angle(p1, a1, p2, a2) < 5.0:
            continue
        pt = geo.ray_intersection(p1, a1, p2, a2)
        if pt is None:
            continue
        xs.append(pt[0])
        ys.append(pt[1])
    if not xs:
        return None
    xs.sort()
    ys.sort()
    return (xs[len(xs) // 2], ys[len(ys) // 2])


def _best_refine_point(valid, area_r, src_est=None, dist_step=100.0,
                       ang_step=10.0, max_dist=1800.0):
    """在区域内搜索补测点 Q，使 Q 处示向度与每条已有示向度的交会角的
    最小值最大（max-min）。

    关键：交会角的评估必须基于**源的位置**。源未知时先用已有示向度估计
    一个 src_est（_bearing_estimate），再假设源就在 src_est，于是 Q 处
    预期示向度 = Q 看向 src_est 的方位。这个代理是自洽的：若估计准确，
    得到的交会角逐个为真；若估计有偏，补测点仍会把三线束的交会角撑开。

    之前用「垂足方向」做代理是错的 —— Q 落在某条示向度线上时方位退化，
    会虚假地给出 90°，实测真实交会角反而只有 1~4°。"""
    if src_est is None:
        src_est = _bearing_estimate(valid)
    if src_est is None:
        return None
    cands = []
    centers = [p for p, _ in valid]
    for cp in centers:
        c = dist_step
        while c <= max_dist:
            a = 0.0
            while a < 360.0:
                q = geo.bearing_point(cp, a, c)
                if math.hypot(q[0], q[1]) <= area_r - 5.0:
                    cands.append(q)
                a += ang_step
            c += dist_step
    best_pt, best_q = None, -1.0
    for q in cands:
        aq = _bearing_to(q, src_est)
        worst = min(geo.crossing_angle(p, pa, q, aq) for (p, pa) in valid)
        if worst > best_q:
            best_q, best_pt = worst, q
            if best_q >= 89.0:
                break
    return best_pt


def adaptive_refine_points(sightings, cleared, area_r=AREA_R,
                           target_angle=80.0, offset_dist=600.0,
                           min_existing_angle=45.0, zero_sig_outer_n=4,
                           offset_ratio=1.4136, offset_dist_min=400.0,
                           offset_dist_max=1700.0):
    """根据粗扫已有的示向度，为交会角不足的频道计算针对性补测点。

    策略：先用已有示向度估计源位置（多交点取中位数），再在区域内搜索一个
    补测点，使其与**每一条**已有示向度的交会角的最小值最大（max-min）。

    为什么不是「垂直偏移 600m」：
      设 P 处测得方位 ang、源距 rho，从 P 沿垂直方向走 c 到 Q，则交会角
      满足 tan(交会角) = c / rho。固定 c=600m 时，rho=1500m 的边缘源只
      得到 atan(600/1500) = 21.8°，恰好压在 estimate_sources 的
      min_angle=20° 门槛上 —— 补测一次（一次路程 + 一次测量）却拿回一个
      仍然病态的解。因此偏离距离必须随 rho 缩放：论文推导的最优交会角
      54.74° 对应 c = rho*tan(54.74°) ≈ 1.414*rho。
      实际实现改为直接对「最小交会角」做网格搜索，比解析式更稳——它能同时
      照顾到「与多条已有示向度都撑开」和「补测点不越界」两个约束。

    只有 1 条示向度的频道（定向源盲区问题）：无法估计源位置，退化为在该
    示向度垂直方向两侧各取一个补测点（距离仍按源距自适应），让至少一个
    能落进定向源的正面半球。

    zero-sightings 频道（粗扫完全没收到信号，定向源朝向背对所有内层点）：
    无法计算垂直方向，改为补 zero_sig_outer_n 个均匀分布的外环点（取 r=area_r-5）。
    这些点以频道组的形式处理：所有 zero-sightings 频道共用同一批外环点，
    每个外环点都测全部 zero-sightings 频道。

    返回三个值：
    - per_channel: {channel: [point, ...]} —— 有示向度频道的专属补测点
    - zero_sig_channels: list[channel] —— 零信号频道列表
    - zero_sig_points: list[point] —— 零信号频道共用的外环补测点
    （调用方需对 zero_sig_points 中每个点测全部 zero_sig_channels）"""
    per_channel = {}
    zero_sig_channels = []

    for c, ss in sightings.items():
        if c in cleared:
            continue
        valid = [(p, a) for (p, a) in ss if a is not None]
        if not valid:
            zero_sig_channels.append(c)
            continue

        if len(valid) == 1:
            # 只有 1 条示向度：垂直方向两侧各取一点，兼顾「撑开交会角」与
            # 「至少一个点落进定向源正面半球」。两侧用同一自适应距离。
            (px, py), ang = valid[0]
            pts = []
            for sign in (+1.0, -1.0):
                rho = _range_along_bearing((px, py), ang, area_r)
                # 注意：变量名不要用 c —— 外层循环的 c 是频道号。
                span = min(max(rho * offset_ratio, offset_dist_min),
                           offset_dist_max)
                qx, qy = geo.bearing_point(
                    (px, py), (ang + sign * 90.0) % 360.0, span)
                d = math.hypot(qx, qy)
                if d > area_r - 5.0:
                    qx, qy = qx * (area_r - 5.0) / d, qy * (area_r - 5.0) / d
                pts.append((qx, qy))
            per_channel[c] = pts
            continue

        best_q = 0.0
        for (p1, a1), (p2, a2) in combinations(valid, 2):
            q = geo.crossing_angle(p1, a1, p2, a2)
            if q > best_q:
                best_q = q

        if best_q >= min_existing_angle:
            continue

        # 已有 >=2 条示向度但交会角不足：先估计源位置，再搜索补测点使
        # 「与每条已有示向度的交会角」的最小值最大（max-min）。
        # 源码曾用「valid[-1] + 固定垂直偏移 600m」，两个问题：
        #   1) 固定距离对边缘源给出的交会角仍 < min_angle（见 docstring）；
        #   2) 补测点裁剪到边界后会进一步压缩夹角。
        best_pt = _best_refine_point(valid, area_r)
        if best_pt is not None:
            per_channel[c] = [best_pt]

    # 零信号频道：外环均匀 zero_sig_outer_n 个点（相位偏移 22.5° 错开内层）
    zero_sig_points = []
    if zero_sig_channels:
        outer_r = area_r - 5.0
        phase_offset = 22.5
        for k in range(zero_sig_outer_n):
            a = math.radians(phase_offset + 360.0 * k / zero_sig_outer_n)
            zero_sig_points.append((outer_r * math.cos(a), outer_r * math.sin(a)))

    return per_channel, zero_sig_channels, zero_sig_points


def two_pass_sweep(client, points, channels=range(1, 21), min_sightings=2,
                   inner_n=7, auto_clear=True, no_signal_limit=None,
                   no_signal_min_r=0.0,
                   early_clear=False, early_clear_min_angle=None,
                   outer_targeted=False, divert=False,
                   divert_min_r=950.0, divert_min_save=150.0):
    """两趟扫描：内层先扫，外层只补内层没测到的频道。

    auto_clear:  扫描时碰到 near 立刻清除，省去后续 homing 整段（默认开启）。
    no_signal_limit: Q3 全向源可传 9，累计无信号达此数就跳过该频道；
        Q4 定向源传 None（无信号可能是背向，不能排除）。
    no_signal_min_r: 仅对距原点 >= 此半径的检测点应用 no_signal_limit。
        Q3 近圈 r=400 的点距边缘源约 1300m（超 1000m 接收半径），那里的
        no_signal 合法，不计入；中远圈 r=1600 才真正反映空频道。
    early_clear: True 时内层扫完立即对已找到的源做 homing 清除，再继续外层
        补盲扫描。省掉先把 18 个外层点全走完再折返清除的往返路程 ——
        内层结束时机器人在区域内部（~1300m 圈），比外层结束时的边缘位置
        （~1795m）更靠近大多数源的估计位置，homing 起始距离平均缩短
        ~300-500m；对 Q3 全向源效果尤著（内层已能发现全部源时外层可全跳）。
    early_clear_min_angle: 仅在 early_clear=True 时生效。传入角度阈值（如
        30.0）时，只对内层已有 ≥2 条示向度且最大交会角 ≥ 此值的频道立即
        homing；其余频道推迟到外层扫完后再清，不打断定向源的补盲扫描。
        传 None（默认）则对所有有示向度的频道立即 homing（Q3 行为）。

    返回 (sightings, cleared)：cleared 是扫描 + early_clear 阶段已清除的频道集合。"""
    cleared = set()
    inner, outer = split_scan_points(points, inner_n)
    origin = getattr(client, "position", (0.0, 0.0)) or (0.0, 0.0)
    inner = order_points_from(origin, inner)
    sightings = detection_sweep(client, inner, channels, min_sightings,
                                auto_clear=auto_clear, cleared=cleared,
                                no_signal_limit=no_signal_limit,
                                no_signal_min_r=no_signal_min_r,
                                divert=divert, divert_min_r=divert_min_r,
                                divert_min_save=divert_min_save)

    if early_clear:
        est = estimate_sources(sightings)
        here = getattr(client, "position", origin) or origin
        if early_clear_min_angle is not None:
            # 有质量门槛：只清内层已有高交会角示向度对的频道，
            # 定向源（内层拿不到足够示向度）自然留给外层补盲。
            eligible = {}
            for c, pos in est.items():
                if c in cleared:
                    continue
                valid = [(p, a) for (p, a) in sightings[c] if a is not None]
                if len(valid) < 2:
                    continue
                best_q = max(geo.crossing_angle(p1, a1, p2, a2)
                             for (p1, a1), (p2, a2) in combinations(valid, 2))
                if best_q >= early_clear_min_angle:
                    eligible[c] = pos
        else:
            # 无质量门槛：对所有有示向度的频道立即 homing（Q3 行为）。
            eligible = {c: v for c, v in est.items() if c not in cleared}
        for ch, x, y in clear_order(here, eligible):
            init_h = [(p, a) for (p, a) in sightings[ch] if a is not None]
            if homing_clear(client, ch, (x, y), init_hist=init_h):
                cleared.add(ch)

    # 外层只测还没集够示向度且未清除的频道
    remain = [c for c in channels
              if c not in cleared and
              sum(1 for s in sightings[c] if s[1] is not None) < min_sightings]
    if remain and outer:
        here = getattr(client, "position", origin) or origin
        if outer_targeted:
            # 用内层已有的示向度估算剩余频道的大致位置，只走能覆盖这些位置的外层点。
            # 0条示向度的频道无法估位（定向源背向全黑），必须保留全部外层点。
            # 只有所有剩余频道都至少有1条示向度时才裁剪。
            seed_sightings = {c: sightings[c] for c in remain}
            zero_signal = [c for c in remain
                           if not any(s[1] is not None for s in sightings[c])]
            if not zero_signal:
                seed_est = estimate_sources(seed_sightings)
                if seed_est:
                    outer_filtered = [p for p in outer
                                      if any(geo.dist(p, pos) <= WORST_RECV
                                             for pos in seed_est.values())]
                    if outer_filtered:
                        outer = outer_filtered
        outer = order_points_from(here, outer)
        extra = detection_sweep(client, outer, remain, min_sightings,
                                auto_clear=auto_clear, cleared=cleared,
                                no_signal_limit=no_signal_limit)
        for c, ss in extra.items():
            sightings[c].extend(ss)
    return sightings, cleared


def estimate_sources(sightings, seed_dist=700.0, min_angle=20.0):
    """对每个检测到信号的频道给出初始位置估计，返回 {channel: (x, y)}。

    只有 1 条示向度的频道不能交会定位（法方程奇异，triangulate 返回 None），
    但绝不能因此丢弃该频道 —— 那等于放弃一个已经发现的干扰源。此时沿这条
    唯一的示向度取 seed_dist 处作为种子点，交给 homing_clear 边走边测去收敛。

    关键：必须按交会角过滤。两条近平行线（交会角 < min_angle）的交点病态，
    偏差可达几千米；用这种解作为 homing 起点，机器人会白白走极远的距离。
    只有找到至少一对交会角合格的示向度时，才用最小二乘定位；否则退回单条
    示向度的射线种子点，让 homing 自己边走边积累高质量交会解。"""
    est = {}
    for c, ss in sightings.items():
        valid = [(p, a) for (p, a) in ss if a is not None]
        pos = None
        if len(valid) >= 2:
            # 找交会角最优的一对（兜底解），同时收集所有质量合格的示向度
            best_score = -1.0
            best_pos = None
            good_set = set()
            for (p1, a1), (p2, a2) in combinations(valid, 2):
                q = geo.crossing_angle(p1, a1, p2, a2)
                if q < min_angle:
                    continue
                p = geo.ray_intersection(p1, a1, p2, a2)
                if p is None or geo.dist(p, (0.0, 0.0)) > 2800:
                    continue
                score = min(geo.dist(p1, p2), 800.0) * math.sin(math.radians(q))
                if score > best_score:
                    best_score = score
                    best_pos = p
                good_set.add((p1, a1))
                good_set.add((p2, a2))
            if best_pos is not None:
                # 只用质量合格的示向度做最小二乘，避免近平行对污染解
                filtered = [(p, a) for (p, a) in valid if (p, a) in good_set]
                ls = geo.triangulate(filtered) if len(filtered) >= 2 else None
                if ls is not None and geo.dist(ls, (0.0, 0.0)) < 2800:
                    pos = ls
                else:
                    pos = best_pos
        if pos is None and valid:
            # 单条示向度或全部近平行：沿最靠近区域中心的检测点取种子
            p, a = min(valid, key=lambda pa: geo.dist(pa[0], (0.0, 0.0)))
            pos = geo.bearing_point(p, a, seed_dist)
        if pos is None and ss:
            pos = ss[0][0]
        if pos is not None:
            est[c] = pos
    return est


def _tour_len(order, origin):
    """巡回总路程：原点 -> 各点 -> 原点。"""
    if not order:
        return 0.0
    L = geo.dist(origin, (order[0][1], order[0][2]))
    for a, b in zip(order, order[1:]):
        L += geo.dist((a[1], a[2]), (b[1], b[2]))
    return L


def two_opt(order, origin, max_pass=40):
    """2-opt 局部搜索改进巡回顺序。

    反复尝试反转任意子段，只要总路程下降就接受，直到无改进或达到 max_pass。
    比近邻贪心稳定省 5%~15% 路程 —— 移动是本题时间的大头，收益直接。"""
    best = list(order)
    best_len = _tour_len(best, origin)
    for _ in range(max_pass):
        improved = False
        for i in range(len(best) - 1):
            for j in range(i + 1, len(best)):
                cand = best[:i] + best[i:j + 1][::-1] + best[j + 1:]
                L = _tour_len(cand, origin)
                if L < best_len - 1e-9:
                    best, best_len, improved = cand, L, True
        if not improved:
            break
    return best


def _scatter(spawn, off, k):
    """以 spawn 为心、黄金角散开的散点，用于失联后重探。"""
    phi = math.radians(137.5 * k)
    return (spawn[0] + off * math.cos(phi), spawn[1] + off * math.sin(phi))


def clear_order(origin, positions, use_two_opt=True):
    """清除顺序：贪心近邻构造 + 起点旋转 + 2-opt + Or-opt(1/2/3) 全程优化。

    复用 order_points_from 保证清除路径与扫描路径用同等质量的优化器。
    对 14-16 个源，Or-opt 能把整块点重新插入到更优位置，比纯 2-opt
    再省约 5-10% 路程。

    返回按顺序排列的 [(channel, x, y), ...]。"""
    items = list(positions.items())
    if not items:
        return []
    pts = [pos for _, pos in items]
    # 位置不唯一时（两个不同源恰好估计在同坐标）保留频道对应关系
    pos_to_ch = {}
    for ch, pos in items:
        pos_to_ch.setdefault(pos, []).append(ch)
    ordered_pts = order_points_from(origin, pts, use_two_opt=use_two_opt)
    result = []
    used = {}
    for pos in ordered_pts:
        chs = pos_to_ch[pos]
        idx = used.get(pos, 0)
        result.append((chs[idx], pos[0], pos[1]))
        used[pos] = idx + 1
    return result


def _mec(points):
    """点集最小外接圆。返回 (圆心, 半径)。

    点数极少（每次交会解 1 个，累计几十个），直接暴搜：最小外接圆必由
    1、2 或 3 个边界点确定，枚举点对的中垂直径圆与三元组外接圆取最小可行者。"""
    pts = [tuple(p) for p in points]
    if not pts:
        return None, 0.0
    if len(pts) == 1:
        return pts[0], 0.0

    def c2(a, b):
        return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2), geo.dist(a, b) / 2

    def c3(a, b, c):
        ax, ay = a; bx, by = b; cx, cy = c
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(d) < 1e-12:
            return c2(*max([(a, b), (b, c), (a, c)], key=lambda pr: geo.dist(*pr)))
        ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay)
              + (cx*cx + cy*cy) * (ay - by)) / d
        uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx)
              + (cx*cx + cy*cy) * (bx - ax)) / d
        return (ux, uy), geo.dist((ux, uy), a)

    def covers(o, r):
        return all(geo.dist(p, o) <= r + 1e-6 for p in pts)

    best = (pts[0], max(geo.dist(pts[0], p) for p in pts))
    for a, b in combinations(pts, 2):
        o, r = c2(a, b)
        if r < best[1] and covers(o, r):
            best = (o, r)
    for a, b, c in combinations(pts, 3):
        o, r = c3(a, b, c)
        if r < best[1] and covers(o, r):
            best = (o, r)
    return best


def homing_clear(client, channel, est, max_probes=48, min_angle=25.0,
                 max_step=250.0, lateral=120.0, clear_try_r=18.0,
                 stop_radius=20.0, init_hist=None):
    """测向—逼近—清除单个源。

    === 为什么不能只沿示向度射线直走 ===
    单条示向度带 ±1° 误差，且误差在空间上连续（同地点重测不变）。实测中
    相隔 250m 的两条示向度交会角仅 0.13°，交会解病态到落在射线外 ——
    沿它走只会振荡，永远进不了 5m 窗口。

    === 本算法 ===
    维护测量历史，每一步：
      1. 当前位置若已 near，就地清除；
      2. 否则在历史中找**交会角 >= min_angle 且基线最长**的一对示向度做交会
         （交会角不足的解一律不采信，这是上一版失败的根因）；
      3. 有可信解就直接走过去；没有就**沿示向度的横向**挪 lateral 拉开基线，
         而不是沿射线往前冲 —— 横向位移才是制造交会角、消除病态的手段。

    init_hist: 扫描阶段已采集的有效示向度列表 [((x,y), svd_deg), ...]，
    直接作为初始 hist 复用，节省开局的横向探测轮次。

    返回是否清除成功。"""
    x, y = est
    # 复用扫描阶段已有示向度作为初始历史，节省开局补测轮次
    hist = list(init_hist) if init_hist else []
    spawn = list(est)            # 起始点，用于限制探索半径
    tried_clear = set()          # 已尝试过清除的位置，避免重复付 3s 失败代价
    lost = 0                     # 连续 no_signal 次数
    last_good = None             # 最近一次有效示向度 (x, y, angle)
    last_sol = None              # 最近一次交会解 —— 失联时优先拿它做清除
    sol_tries = 0                # 围绕交会解的微偏移重试计数
    since_sol = 99               # 距上次算出交会解经过了几次测量
    micro = 20.0                 # 近距离失联后的微小前进量
    visits = {}                  # 位置 -> 到访次数，用于识别极限环
    sols = []                    # 历次交会解，用于最小外接圆停机判据
    # 失联后的恢复偏移：散点重探。定向源只有 180° 覆盖锥，单个点落在锥外
    # 返回的 no_signal 与「无源」「超距」不可区分，只能换位置重探。
    # 最后的 None 是「近距离失联」专用分支（见主循环）。
    recover = (180.0, 130.0, 250.0, 360.0, 450.0, 540.0, 700.0, None, None, None)

    def _try_clear(px, py):
        """在 (px,py) 处尝试清除；同一位置只试一次。"""
        key = (round(px, 0), round(py, 0))
        if key in tried_clear:
            return False
        tried_clear.add(key)
        return client.clear(px, py, channel).get("clear_result") == "success"

    # 初始估计来自多条示向度的最小二乘，是清除半径内最可信的位置之一。
    # 先拿它试一次：定向源的锥外盲区会让这一步返回 no_signal（测不到方向），
    # 但**清除只看光学对准、与覆盖锥无关**。不先试这一下，就会在没有任何
    # 交会解的情况下散点跳走，白白错过一个已经在 20m 内的解。
    if _try_clear(x, y):
        return True
    # 不把初始估计加入 sols —— 它可能是近平行示向度的病态解，
    # 污染 MEC 停机判据会让 homing 过早或过晚终止。
    # sols 只收录 homing 过程中实际测量得到的交会解。

    for _ in range(max_probes):
        # 极限环检测：位置被反复到访说明状态机在原地打转（实测过三位置
        # A->B->C->A 的死循环，每次都从「沿横向挪 lateral」回到同一点）。
        # 一旦发现，就把这个位置拉黑并强制走出去，否则会跑满 max_probes。
        k = (round(x, 0), round(y, 0))
        visits[k] = visits.get(k, 0) + 1
        if visits[k] > 2:
            ang = math.radians(137.5 * visits[k])
            r_esc = 250.0 + 90.0 * visits[k]
            x, y = _scatter(spawn, r_esc, visits[k])
            if _try_clear(x, y):
                return True
            continue
        r = client.measure(x, y, channel)
        mr = r.get("measure_result")
        if mr == "near":
            return _try_clear(x, y)
        if mr != "direction":
            # no_signal：锥外 / 超距 / 无源，三者不可区分。
            # 判据：如果刚拿到过示向度就走到了这里，说明这一步迈到了源附近
            # 或锥边缘 —— 优先试清除，而不是散点跳走（跳走会退回到已测点，
            # 造成两点之间的死循环振荡，实测每次跑满 max_probes）。
            if last_sol is not None and since_sol <= 1:
                off = None
            elif lost >= len(recover):
                return False
            else:
                off = recover[lost]
            if off is None:
                # 「近距离失联」分支：已经走到源上（几米内）却因锥朝向拿到
                # no_signal。此时退回远处重测毫无意义——交会解比任何单条
                # 示向度都准，应当沿「最近一次有效示向度」方向小步推进，
                # 逐次尝试清除，直到走出这一步的合理范围。
                # (1) 先试最近一次交会解：它比任何单条示向度都准。定向源在
                #     锥外会返回 no_signal，但「走到解上再清除」是有效的 ——
                #     清除要求的是光学对准，与覆盖锥无关。
                if last_sol is not None and sol_tries < 8:
                    # 在交会解周围做微小偏移重试：解本身通常已落在源附近
                    # （实测 0.3m），1~2m 的偏移足以跨过 20m 清除边界。
                    sx, sy = last_sol
                    off = sol_tries * 2.0
                    phi = math.radians(137.5 * sol_tries)
                    sol_tries += 1
                    if _try_clear(sx + off * math.cos(phi),
                                  sy + off * math.sin(phi)):
                        return True
                    continue
                # (2) 没有交会解：沿最近一次有效示向度小步推进再试
                if last_good is not None and micro <= 160.0:
                    gx, gy, ga = last_good
                    x, y = geo.bearing_point((gx, gy), ga, micro)
                    micro += 35.0
                    if _try_clear(x, y):
                        return True
                    continue
                return False
            lost += 1
            # 有交会解时从解出发散点，而非从初始估计点 spawn：
            # 初始估计可能因近平行示向度而偏几百米，固定从 spawn 散点
            # 会永远打不到真实位置附近。
            center = last_sol if last_sol is not None else tuple(spawn)
            x, y = _scatter(center, off, lost)
            continue
        lost = 0
        hist.append(((x, y), r["svd_deg"]))
        a = r["svd_deg"]
        last_good = (x, y, a)
        since_sol += 1

        # 选交会角合格的一对。按「交会角的定位精度」而非「基线长度」排序：
        # 交会解的水平误差 ~ base_eff * sin(dθ) / sin(gamma)，真正决定精度的是
        # 交会角 gamma，基线只在 gamma 相当时才起作用。此前按最长基线选，
        # 实测选中过一个交会角 38°/基线 1218m、距真值 18.7m 的解，而同一批
        # 数据里存在交会角 73.38°/基线 570m、距真值 1.1m 的解被 base>=40 的
        # 门槛丢掉 —— 那一对的两条线来自相隔 16.7m 的两点，正是上一步刚走出来的。
        #
        # 近距离放宽角度门槛：机器已在源附近（有交会解且距解 < 200m）时，
        # ±1° 误差下交会角 15° 仍能给出 < 70m 残差，足以落入清除半径。
        # 远距离保持 min_angle 以避免病态解。
        if last_sol is not None and geo.dist((x, y), last_sol) < 200.0:
            eff_min_angle = max(15.0, min_angle * 0.6)
        else:
            eff_min_angle = min_angle
        best = None
        for (px, py), pa in hist[:-1]:
            base = geo.dist((px, py), (x, y))
            q = geo.crossing_angle((px, py), pa, (x, y), a)
            if q < eff_min_angle:
                continue
            score = min(base, 600.0) * math.sin(math.radians(q))
            if best is None or score > best[0]:
                p = geo.ray_intersection((px, py), pa, (x, y), a)
                if p is not None:
                    best = (score, p)

        if best is not None:
            tx, ty = best[1]
            last_sol = (tx, ty)
            since_sol = 0
            sols.append((tx, ty))

            # 停机判据（几何，而非启发式）：历次交会解的**最小外接圆**半径
            # 就是当前对源位置的不确定度。≤ 清除半径 20m 时，圆内任意点都在
            # 某个源的清除范围内，此时清除的期望收益已经饱和 —— 再补测只会
            # 白付 5s/次的测量成本。反过来，半径 > 20m 时就算单次交会解看着
            # 很近也不可靠（±1° 误差下解本身可能有 5~15m 残差）。
            # 取圆心而非最近一次解：圆心是全部证据的折中，比任何单条示波稳定。
            if len(sols) >= 2:
                co, cr = _mec(sols)
                if co is not None and cr <= stop_radius:
                    if _try_clear(co[0], co[1]):
                        return True
                    # 圆心清除失败：说明解整体偏了，继续补测收缩不确定度
            d = geo.dist((x, y), (tx, ty))
            if d <= clear_try_r:
                # 已进入清除半径：直接尝试清除。交会解受 ±1° 误差影响必然
                # 有 5~15m 残差，非要等到 5m 的 near 窗口只会白跑冤枉路。
                if _try_clear(tx, ty):
                    return True
                # 清除失败：说明解偏了，改为横向拉基线继续测
                x, y = geo.bearing_point((x, y), (a + 90.0) % 360, lateral)
                continue
            step = min(d, max_step)
            x, y = geo.bearing_point((x, y),
                                     math.degrees(math.atan2(ty - y, tx - x))
                                     % 360, step)
            continue

        # 无可信交会解：沿横向挪一步拉开基线（这一步是收敛的关键）
        # 自适应步长：距源估计越远，横向步需要越大才能撑开足够交会角。
        # 有交会解时用解距离估算，否则用到初始估计点的距离兜底。
        if last_sol is not None:
            est_dist = geo.dist((x, y), last_sol)
        else:
            est_dist = geo.dist((x, y), tuple(spawn))
        # 目标交会角 30°：lateral = dist * tan(30°) ≈ dist * 0.577，
        # 下限 80m（近距离不要步子太大），上限 350m（避免跑出接收范围）。
        adaptive_lateral = max(80.0, min(350.0, est_dist * 0.577))
        x, y = geo.bearing_point((x, y), (a + 90.0) % 360, adaptive_lateral)

    return False
