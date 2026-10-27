"""
问题 3 / 问题 4 主程序：干扰源的自动定位与清除。

用法：
1. 先启动模拟器并开始一次测试，等待界面显示机器狗接口就绪。
2. python robot.py                      # Q3（全向源），使用演示标识
   python robot.py --q4                 # Q4（含定向源），7 点内层 + 自适应精扫
python robot.py --id <YOUR_TEAM_ID>  # 真实模拟器使用自己的队号
   python robot.py --id test            # 本地 mock 联调

本地 A/B 用的环境变量（只在 robot.py 里读，模拟器侧无影响）：
   CUMCM_DIVERT=1   打开顺路分流（默认关，实测净负，仅用于复现实验）
"""
import argparse
import os
import time

import strategy
from client import RobotClient, NoResponse

DEFAULT_ROBOT_ID = os.environ.get("ROBOT_ID", "test")


def main():
    parser = argparse.ArgumentParser(description="CUMCM 2026 B 题机器狗主程序")
    parser.add_argument("--id", default=DEFAULT_ROBOT_ID,
                        help="模拟器标识（默认 %(default)s，可用 ROBOT_ID 覆盖）")
    parser.add_argument("--port", type=int, default=2026,
                        help="模拟器端口（默认 %(default)s）")
    parser.add_argument("--q4", action="store_true",
                        help="Q4 模式：7 点内层 + 自适应精扫，不启用 no_signal_limit")
    args = parser.parse_args()

    # 逐动作轨迹落盘（官方加密日志无法解读，复盘只能靠这个）。
    # 每局一个文件，文件名带模式与时间戳，便于按局对比。
    _trace = "trace_%s_%s.jsonl" % (
        "q4" if args.q4 else "q3",
        time.strftime("%m%d_%H%M%S"))
    client = RobotClient(args.id, base_url=f"http://127.0.0.1:{args.port}",
                         trace_path=_trace)
    print("轨迹文件:", _trace)

    # 1. 进入（倒计时/接口未开放期间连接会失败，轮询重试）
    r = None
    for _ in range(600):
        try:
            r = client.enter()
            break
        except NoResponse:
            time.sleep(1.0)
    if r is None:
        print("无法进入目标区域")
        return
    print("进入成功，可用现实时间(s):", r.get("remaining_real_duration_s"))

    # 2. 覆盖扫描
    #    Q3：7 点（圆心 + 均匀环 r=1122.96），no_signal_limit=18 快速跳过空频道
    #    Q4：7 点（圆心 + 均匀环 r=1122.96），no_signal_limit=None
    inner_n = 7
    no_sig_min_r = 0.0
    if args.q4:
        # Q4 定向源：三阶段自适应精扫。
        # 阶段1：7点内层粗扫，收集每个频道的初始方位角。
        # 阶段2：adaptive_refine_points 对交会角 <45° 的频道补专属垂直点。
        # 阶段3：全量示向度汇总后 Or-opt 顺序统一清除。
        # 注：零信号频道（定向源朝向背对所有内层点）无法廉价补测，
        #     实验表明加 4 个外环点路程增加 4000m+ 收益不足，保持 7 点策略。
        #
        # inner_r=1122.96（= Q3 用值），显式传以防 coverage_points 默认值回退。
        # 曾漏传此参数，用了当时的默认值 area_r - worst_recv/2 = 1300：
        #   r=1300  相邻检测点间距 1300m，环周长 7800m，最坏内部覆盖 936.5m
        #   r=1123  相邻检测点间距 1123m，环周长 6738m，最坏内部覆盖 1000.0m
        # 1123 更短（每步省 177m，6 步省 1062m）且覆盖恰好卡在接收半径下界。
        # 实测 Q4 三局都在走 1300m 的环步（1300/1300/1300/1300/1300），
        # 前 5 个点常常拿不到信号也照走 —— 每步白多付 177m。
        # （coverage_points 的默认值现已改成同一个闭式最优解，见其 docstring。）
        points = strategy.coverage_points(inner_r=1122.96)[:inner_n]
        no_sig_limit = None
        do_early_clear = False
        early_clear_angle = None
        do_outer_targeted = False
        do_divert = False
    else:
        # Q3 全向源：内层 = 圆心 + 6 点均匀环 r=1122.96。
        # 该半径是「7 点按最坏接收半径 1000m 覆盖全区」的**最小可行值**
        # （最坏覆盖随 r 单调递减，1122.96 处恰好等于 1000m），因此也是最
        # 省路程的环：周长 6738m。推导见 strategy.coverage_points docstring。
        # 旧的「400×3 + 1600×4」混合布局已废弃 —— 它在 (1270,-1270) 一带有
        # 1312m 的覆盖洞，且边缘源在 recv=1250 时退化为单点无法定位。
        #
        # === 只用内层 7 点，不加外环（2026-09-11 实测后定案）===
        # 曾试过内层 7 点 + 外环 6 点（r=1600，相位错开半格）。纯几何扫描支持
        # 它：13 点巡回 10385m vs 7 点 6738m，但能定位的源比例从 92.7% 升到
        # 98.0%。**实测否掉了它**：那一局外层一个频道都没需要补测（内层 7 点
        # 已发现全部源），13 点却多走了 3647m = 729s 的确定开销，vt/源 从
        # 338.0 涨到 451.2（+33.5%）。
        #
        # 根因：全向源在内层 r=1122.96 已被 100% 覆盖（区域内最大最近距离恰好
        # 1000m，是 Kershner 下界的临界值），外环对全向源没有任何发现价值。
        # 外环唯一能救的是「定向源朝区域外、内层 7 点全落在 180° 盲锥里」，
        # 而官方附件说明**全向/定向是每个源的属性（Q3/Q4 都含两类）**，
        # 朝向随机 —— 实测多局未出现这种源。
        #
        # 若将来漏源（检测数 < 实际源数）且确认是定向源，再加外环：
        #   strategy.coverage_points(inner_r=1122.96, outer_r=1600.0,
        #                            outer_n=6, outer_phase=30.0)
        # 并把 inner_n 改错会让外环失效：inner_n 必须是**内层点数 7**，
        # 它同时控制 points[:inner_n] 切片和 split_scan_points 的内/外分层。
        #
        # early_clear_min_angle=30.0：只有交会角达标的源才在内层扫完后立即清，
        # 避免角度差的估计触发 homing 白跑（之前 None 导致 4 次失败）。
        # 注意切片：coverage_points 默认还带 outer_n=18 的外环，必须取前 7 点。
        points = strategy.coverage_points(inner_r=1122.96)[:inner_n]
        # 空频道跳过限值：内层环 r=1123 到边缘源约 700~1000m，基本都在接收
        # 半径内，连续无信号才真正指向空频道。降到 4/9 会误伤边缘源，保持宽松。
        no_sig_limit = 18
        no_sig_min_r = 0.0
        do_early_clear = True
        early_clear_angle = 30.0
        do_outer_targeted = False
        # 顺路分流：**默认关闭**，见 strategy.divert_worth 的说明。
        # 本地 A/B 用 CUMCM_DIVERT=1 打开（不改代码）。
        #
        # 为什么不默认开（2026-09-12 同种子 A/B 实测）：
        #   种子 42  -285s(-6.1%)   种子 21  -265s(-5.4%)   种子 13  -13s(-0.3%)
        #   种子  7   +59s(+1.3%)   种子 99 +1271s(+24.1%)   <- 严重回退
        # 均值净负。根因不是判据松，是判据**算错了对象**：divert_worth 把
        # 「从 p_now 走到 p_wait」这块**必须走**的环上路程当成了「省下的」，
        # 于是源越远、越在岛外，判据越倾向于分流 —— 方向正好反了。
        # 真正的代价是「绕出去再绕回来」的纯绕行差，实测 ch20=271m ch6=337m
        # （种子 99 那局，在检测点 ±150m 噪声下完全无法区分），而分流一次就把
        # 内层环的连续遍历打断：OFF 走的是干净环 1123m×6，ON 变成
        # 1686/2337/2274/2937/2570 的折返段，环遍历本身 6738m -> 14413m。
        do_divert = os.environ.get("CUMCM_DIVERT", "0") == "1"
    print(f"覆盖检测点: {len(points)}（{'Q4' if args.q4 else 'Q3'}）")
    sightings, cleared = strategy.two_pass_sweep(
        client, points, auto_clear=False, min_sightings=2,
        no_signal_limit=no_sig_limit, no_signal_min_r=no_sig_min_r,
        inner_n=inner_n,
        early_clear=do_early_clear,
        early_clear_min_angle=early_clear_angle,
        outer_targeted=do_outer_targeted,
        divert=do_divert)
    print("粗扫阶段已清除:", sorted(cleared))
    found = [c for c, ss in sightings.items() if ss]
    print("粗扫检测到信号的频道:", sorted(found))

    # 2.5 Q4 自适应精扫：对交会角不足的频道计算垂直补测点，只走这些点
    if args.q4:
        # zero_sig_outer_n：零信号频道的外环补测点数。
        #
        # === 为什么从 0 改成 6 ===
        # 零信号频道在 Q4 里**不是空频道，是漏掉的定向源**。本地 mock 带真值
        # 统计（_att/miss.py，4 个种子 × 14 源 × 50% 定向）：
        #     seed 5 漏 1 个 ch13(定向 r=1335 从未收到信号)
        #     seed 7 漏 2 个 ch6(r=1668) ch13(r=1607)
        #     seed 11 漏 2 个 ch19(r=1697) ch20(r=1749)
        #     seed 23 漏 3 个 ch1(r=1738) ch2(r=1323) ch8(r=1688)
        # 漏源率 7~21%（均值 ~14%），**8/8 全是定向源、全 r>1300、
        # 全部一条信号都没收到**。
        #
        # 几何上这是必然的：定向源的盲区是它**背面**的半平面，源在 r≈1700
        # 且朝向区域外时，可见半圆与 r<=1800 圆盘的交集整个落在内层环
        # r=1123 之外 —— 7 个内层点没有一个落在里面，无论怎么调内层参数。
        # 而 adaptive_refine_points 对 0 条示向度的频道只能进 zero_sig 分支，
        # zero_sig_outer_n=0 时**一个补测点都不生成**，等于直接放弃。
        #
        # 外环补测点可见率（源 r∈[1400,1795]、朝向朝外±15°，点铺在 r=1795）：
        #     n=0 ->   0%    n=4 -> 55.7%   n=6 -> 81.6%   n=8 -> 91.8%
        #     n=12 -> 97.6%  n=18 -> 99.6%
        # 取 3（不是 6）：本题按「时间为主，漏源只作为论文里的缺陷记录」计分。
        # 7 种子 × (14 源, 50% 定向) A/B 实测（_att/sweep_outer.py）：
        #     n=0  漏源 1.71/12   vt 6057/6186
        #     n=2  漏源 1.00/7    vt 6929/7002
        #     n=3  漏源 1.00/7    vt 6771/6705   <- 拐点
        #     n=4  漏源 0.57/4    vt 7387/7616
        #     n=6  漏源 0.20/1    vt 8020/7982
        #     n=8  漏源 0.20/1    vt 8085/8233
        # 3 -> 4 要多付 616s 才少漏 0.43 个源；3 -> 6 多付 1249s 才少漏 0.8 个。
        # 时间为主的前提下，拐点 3 是最优点。
        #
        # 注意这**只在 zero_sig_channels 非空时才生成**，全向源局（Q3）
        # 和没有零信号频道的 Q4 局完全不付这个代价。
        # 环境变量 CUMCM_OUTERN 覆盖外环点数（默认 3），用于本地 A/B。
        # Q3 全向源：零信号频道就是空频道，外环补测无意义，置 0 跳过。
        # Q4 定向源：零信号可能是背向内层的定向源，需要外环探测。
        if args.q4:
            zero_outer_n = int(os.environ.get("CUMCM_OUTERN", "3"))
        else:
            zero_outer_n = 0
        per_channel, zero_sig_channels, zero_sig_points = \
            strategy.adaptive_refine_points(sightings, cleared,
                                            zero_sig_outer_n=zero_outer_n)
        # 建立 point -> 需要在该点测的频道列表
        # 有示向度的频道：每个点只测对应的特定频道
        # 零信号频道：在每个外环补测点上测全部 zero_sig_channels
        pt_to_channels = {}
        all_refine_pts = []
        for c, pts in per_channel.items():
            for pt in pts:
                if pt not in pt_to_channels:
                    pt_to_channels[pt] = []
                    all_refine_pts.append(pt)
                pt_to_channels[pt].append(c)
        for pt in zero_sig_points:
            if pt not in pt_to_channels:
                pt_to_channels[pt] = []
                all_refine_pts.append(pt)
            pt_to_channels[pt].extend(zero_sig_channels)
        if zero_sig_channels:
            print(f"零信号频道 {sorted(zero_sig_channels)}：补 {len(zero_sig_points)} 个外环点")
        if all_refine_pts:
            here = client.position or (0.0, 0.0)
            ordered_refine = strategy.order_points_from(here, all_refine_pts)
            print(f"精扫补测点: {len(ordered_refine)} 个，涉及频道: "
                  f"{sorted(set(c for cs in pt_to_channels.values() for c in cs))}")
            for pt in ordered_refine:
                for c in pt_to_channels.get(pt, []):
                    if c in cleared:
                        continue
                    # 已经够 min_sightings 条有效示向度就不必再测。
                    # detection_sweep 里有同样的跳过逻辑，但精扫这段一直没有，
                    # 于是零信号频道在外环点上被反复测：实测一局 ch1 被测 11
                    # 次（第 4 次就拿到示向度了，之后又白测 7 次 ×5s）。
                    # 外环点每个点要测全部零信号频道，这个浪费是乘上去的。
                    if sum(1 for s in sightings.get(c, []) if s[1] is not None) \
                            >= 2:
                        continue
                    r = client.measure(pt[0], pt[1], c)
                    mr = r.get("measure_result")
                    if mr == "direction":
                        sightings[c].append((pt, r["svd_deg"]))
                    elif mr == "near":
                        sightings[c].append((pt, None))
                        if client.clear(pt[0], pt[1], c).get("clear_result") == "success":
                            cleared.add(c)
        else:
            print("精扫：所有频道交会角已足够，无需补测")

    # 3. 定位
    est = strategy.estimate_sources(sightings)
    for c, (x, y) in sorted(est.items()):
        if c in cleared:
            continue
        print(f"  频道{c}: ({x:.1f}, {y:.1f})")

    # 4. 近邻顺序清除（跳过扫描阶段已清除的频道）
    here = client.position or (0.0, 0.0)
    order = strategy.clear_order(here, {c: v for c, v in est.items() if c not in cleared})
    cleared_count = len(cleared)
    for ch, x, y in order:
        init_hist = [(p, a) for (p, a) in sightings.get(ch, []) if a is not None]
        # 标记阶段：事后可从轨迹里切出每个频道的 homing 全过程，
        # 看 fail 发生在哪个分支（交会解偏差 / 走完 recover 全序列）。
        client.note(f"homing:ch{ch}:est=({x:.0f},{y:.0f}):hist={len(init_hist)}")
        ok = strategy.homing_clear(client, ch, (x, y), init_hist=init_hist)
        if ok:
            cleared_count += 1
            print(f"清除 频道{ch} 成功")
        else:
            print(f"清除 频道{ch} 失败（需复盘）")

    # 5. 退出
    client.note("exit")
    r = client.exit()
    vt = round(client.virtual_time, 2)
    print("退出:", r.get("exit_reason"), "| 虚拟时间:", vt,
          "| 清除:", cleared_count, "个")
    if cleared_count:
        print("平均定位清除时间 =", round(vt / cleared_count, 2), "秒")
    # 时间分项分析
    t_move  = round(client.total_dist / 5, 1)
    t_meas  = round(client.n_measure * 5, 1)
    t_sw    = round(client.n_channel_switch * 1, 1)
    t_clear = round(client.n_clear_ok * 3 + client.n_clear_fail * 5, 1)
    print(f"  移动: {t_move}s ({round(t_move/vt*100)}%)  "
          f"测量×{client.n_measure}: {t_meas}s ({round(t_meas/vt*100)}%)  "
          f"切换×{client.n_channel_switch}: {t_sw}s  "
          f"清除ok={client.n_clear_ok}/fail={client.n_clear_fail}: {t_clear}s  "
          f"总路程: {round(client.total_dist)}m")
    client.close()
    print("轨迹已保存:", _trace)


if __name__ == "__main__":
    main()
