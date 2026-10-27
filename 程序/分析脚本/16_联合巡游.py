# -*- coding: utf-8 -*-
"""把「精扫补测点」和「清除目标」合成一条路径的排程器。

动机（Q4 实测，trace_q4_0912_024117）：
    内层 7 点结束 vt=2025，之后精扫段走到 vt=6833（+4808s，占全 66%），
    其中精扫段路程 19185m（=3837s），最后清除尾段又要 770m。
    两段是**分开排的**：精扫段走完 36 个补测点后，清除阶段再从 robots 停下的
    位置重新贪心排序去清 14 个源。但实际上精扫点里有相当一部分就在源附近
    —— 扫到那儿顺手清掉，比走完全部再折返便宜。

本模块只做排程，不改测量协议：给一组「点节点（去测某个频道）」和一组
「清除节点（去清某个频道）」，返回一条统一的访问顺序，相邻节点共用一次移动。

用法：
    seq = plan(origin, measure_nodes, clear_nodes)
    # seq: [("measure", pt, channel), ("clear", pt, channel), ...]
"""
import math
import sys

sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from strategy import geo  # noqa: E402


def _nn(origin, pts, key=None):
    """最近邻顺序。"""
    rem = list(pts)
    out = []
    cur = origin
    while rem:
        i = min(range(len(rem)), key=lambda k: geo.dist(cur, rem[k]))
        out.append(rem.pop(i))
        cur = out[-1]
    return out


def _tour_len(order, origin):
    if not order:
        return 0.0
    L = geo.dist(origin, order[0])
    for i in range(1, len(order)):
        L += geo.dist(order[i - 1], order[i])
    return L


def _two_opt(order, origin, max_pass=40):
    if len(order) < 3:
        return order
    n = len(order)
    best = list(order)
    best_len = _tour_len(best, origin)
    improved = True
    passes = 0
    while improved and passes < max_pass:
        improved = False
        passes += 1
        for i in range(n - 1):
            for j in range(i + 2, n):
                cand = best[:i + 1] + best[i + 1:j + 1][::-1] + best[j + 1:]
                L = _tour_len(cand, origin)
                if L < best_len - 1e-9:
                    best, best_len = cand, L
                    improved = True
    return best


def plan(origin, measure_nodes, clear_nodes, merge_dist=420.0):
    """返回合并后的执行顺序。

    measure_nodes: [(pt, channel), ...]  —— 精扫补测点
    clear_nodes:   [(pt, channel), ...]  —— 估计出的源位置

    合并判据：若某个清除节点与某个补测点的距离 <= merge_dist，就把清除
    插到该补测点**之后**执行 —— 机器人反正要走到那里，顺手清掉，省下
    后来从场另一头折返的整段穿场。

    未被合并的清除节点仍在末尾按最近邻排序（这些源离补测点太远，先扫完
    再统一清更省）。

    返回 [("measure", pt, ch) | ("clear", pt, ch), ...]。clear 节点的 pt 用
    源估计位置（homing 会边走边修正）。
    """
    if not measure_nodes:
        order = _nn(origin, [pt for pt, _ in clear_nodes])
        bypt = {pt: ch for pt, ch in clear_nodes}
        return [("clear", pt, bypt[pt]) for pt in order]

    mpts = [pt for pt, _ in measure_nodes]
    morder = _nn(origin, mpts)
    morder = _two_opt(morder, origin)

    # 补测点 -> 频道
    pt2ch = {}
    for pt, ch in measure_nodes:
        pt2ch.setdefault(pt, []).append(ch)

    seq = [("measure", pt, None) for pt in morder]   # channel 在执行时按 pt2ch 取
    joined = []
    leftover = []
    for cpt, cch in clear_nodes:
        # 找最近的补测点
        best_i, best_d = None, None
        for i, mpt in enumerate(morder):
            d = geo.dist(cpt, mpt)
            if best_d is None or d < best_d:
                best_i, best_d = i, d
        if best_d is not None and best_d <= merge_dist:
            joined.append((best_i, best_d, cpt, cch))
        else:
            leftover.append((cpt, cch))

    # 按补测点序号插入（同一补测点后可能挂多个清除，按距离近的先清）
    joined.sort(key=lambda t: (t[0], t[1]))
    out = []
    ins = {}
    for i, d, cpt, cch in joined:
        ins.setdefault(i, []).append((cpt, cch))
    for i, item in enumerate(seq):
        out.append(item)
        for cpt, cch in ins.get(i, []):
            out.append(("clear", cpt, cch))

    lorder = _nn(out[-1][1] if out else origin, [pt for pt, _ in leftover])
    lmap = {pt: ch for pt, ch in leftover}
    for pt in lorder:
        out.append(("clear", pt, lmap[pt]))

    # 把 measure 节点的 channel 补上（同一个点可能测多个频道，这里每个点
    # 只放一个占位，执行端再展开）
    return out, pt2ch
