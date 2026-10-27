"""
多案例批量测试：对多个随机种子跑完整策略，统计清除比例与平均定位清除时间。

每个种子起一个独立 mock 进程（独立端口、独立状态），跑完打印汇总表。
用法：
    python batch_test.py                 # 默认 8 个种子，全向
    python batch_test.py --directional 0.4   # 40% 定向源（问题 4 预演）
"""
import argparse
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

import strategy
from client import RobotClient, NoResponse


def wait_ready(port, tries=60):
    """轮询直到 mock 接口可用。"""
    url = f"http://127.0.0.1:{port}/enter"
    for _ in range(tries):
        try:
            urllib.request.urlopen(
                urllib.request.Request(
                    url, data=b"{}", method="POST",
                    headers={"Content-Type": "application/json"}),
                timeout=1.0)
            return True
        except urllib.error.HTTPError:
            return True          # 有响应即说明端口已开放
        except Exception:
            time.sleep(0.3)
    return False


def parse_truth(text):
    """从 mock stdout 解析真实干扰源，返回 {channel: (x, y, 是否定向)}。"""
    truth = {}
    for m in re.finditer(
            r"频道(\d+): \(([-\d.]+), ([-\d.]+)\).*?(全向|定向)", text):
        truth[int(m.group(1))] = (float(m.group(2)), float(m.group(3)),
                                  m.group(4) == "定向")
    return truth


def run_case(seed, port, n_sources, directional):
    proc = subprocess.Popen(
        [sys.executable, "-u", "mock_simulator.py",
         "--seed", str(seed), "--port", str(port),
         "--sources", str(n_sources), "--directional", str(directional)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace")
    try:
        if not wait_ready(port):
            return None
        client = RobotClient("test", base_url=f"http://127.0.0.1:{port}")
        client.enter()

        pts = strategy.coverage_points()
        sightings, _ = strategy.two_pass_sweep(client, pts)
        est = strategy.estimate_sources(sightings)
        order = strategy.clear_order((0.0, 0.0), est)

        cleared = 0
        failed = []
        for ch, x, y in order:
            if strategy.homing_clear(client, ch, (x, y)):
                cleared += 1
            else:
                failed.append(ch)
        vt = client.virtual_time
        client.exit()
        return {
            "seed": seed, "found": len(est), "cleared": cleared,
            "failed": failed, "vt": vt,
            "avg": vt / cleared if cleared else float("inf"),
            "n_actions": len(client.log),
        }
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--sources", type=int, default=13)
    ap.add_argument("--directional", type=float, default=0.0)
    ap.add_argument("--port0", type=int, default=2100)
    args = ap.parse_args()

    rows = []
    for i in range(args.seeds):
        r = run_case(i + 1, args.port0 + i, args.sources, args.directional)
        if r is None:
            print(f"seed {i+1}: mock 启动失败")
            continue
        rows.append(r)
        print(f"seed {r['seed']:>2}  发现 {r['found']:>2}/{args.sources}  "
              f"清除 {r['cleared']:>2}/{args.sources}  "
              f"平均 {r['avg']:>8.1f}s  动作 {r['n_actions']:>4}  "
              f"失败频道 {r['failed']}")

    if rows:
        n = len(rows)
        tot_c = sum(r["cleared"] for r in rows)
        print("\n==== 汇总 ====")
        print(f"案例数 {n}")
        print(f"清除比例 {tot_c}/{n*args.sources} = "
              f"{tot_c/(n*args.sources)*100:.1f}%")
        finite = [r["avg"] for r in rows if r["avg"] != float("inf")]
        if finite:
            print(f"平均定位清除时间 {sum(finite)/len(finite):.1f}s  "
                  f"(最好 {min(finite):.1f}s / 最差 {max(finite):.1f}s)")


if __name__ == "__main__":
    main()
