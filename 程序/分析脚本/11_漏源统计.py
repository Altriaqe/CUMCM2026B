# -*- coding: utf-8 -*-
"""在已知真值的 mock 世界上统计 Q4 的**漏源率**。

官方模拟器不给真值，所以「漏了几个源」这件事在正式测试里看不出来 ——
只能靠本地 mock（带 /_truth 端点）统计。这是唯一能验证「零信号频道」
到底是空频道还是漏源的途径。

用法：
    python _att/miss.py <port> <seed> <sources> <directional> [q3|q4]

脚本自己起停 mock，每局结束后拉 /_truth 与 robot.py 的输出对比：
    - 真实源数 / 定向源数
    - robot 检测到的频道数
    - robot 清除成功的频道数
    - **漏源**：真实存在但没被清除的频道（分类：从未收到信号 / 收到了但没清）
"""
import json
import math
import os
import re
import socket
import subprocess
import sys
import urllib.request

sys.stdout = __import__('io').TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable


def wait_port(port, timeout=20.0):
    import time
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            s = socket.create_connection(("127.0.0.1", port), timeout=0.4)
            s.close()
            return True
        except OSError:
            time.sleep(0.3)
    return False


def kill_port(port):
    try:
        out = subprocess.run(["netstat", "-ano"], capture_output=True, text=True,
                             timeout=20).stdout
    except Exception:
        return
    for line in out.splitlines():
        if "127.0.0.1:%d " % port in line and "LISTENING" in line:
            subprocess.run(["taskkill", "/F", "/PID", line.split()[-1]],
                           capture_output=True)


def truth(port):
    try:
        d = json.loads(urllib.request.urlopen(
            "http://127.0.0.1:%d/_truth" % port, timeout=5).read())
    except Exception:
        return None
    return d.get("sources") if isinstance(d, dict) else d


def one(port, seed, nsrc, direc, mode, outer=None):
    kill_port(port)
    import time
    time.sleep(0.5)
    mock = subprocess.Popen(
        [PY, "-u", "mock_simulator.py", "--port", str(port), "--seed", str(seed),
         "--sources", str(nsrc), "--directional", str(direc)],
        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not wait_port(port):
        print("mock 起不来")
        return None
    try:
        cmd = [PY, "-u", "robot.py", "--id", "miss", "--port", str(port)]
        if mode == "q4":
            cmd.append("--q4")
        env = dict(os.environ)
        if outer is not None:
            env["CUMCM_OUTERN"] = outer
        pr = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                            timeout=3600, env=env)
        out = pr.stdout
        t = truth(port)
    finally:
        mock.terminate()
        try:
            mock.wait(timeout=5)
        except Exception:
            mock.kill()

    real = {s["channel"]: s for s in (t or [])}
    # robot 侧：从 stdout 抓「粗扫检测到信号的频道」与「清除 频道N 成功」
    found = set()
    m = re.search(r"粗扫检测到信号的频道: \[([^\]]*)\]", out)
    if m:
        found = {int(x) for x in m.group(1).split(",") if x.strip()}
    cleared = {int(x) for x in re.findall(r"清除 频道(\d+) 成功", out)}
    vt = re.search(r"虚拟时间: ([\d.]+)", out)
    return {
        "real": real, "found": found, "cleared": cleared,
        "vt": float(vt.group(1)) if vt else None,
        "raw": out,
    }


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 2060
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    nsrc = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    direc = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    mode = sys.argv[5] if len(sys.argv) > 5 else "q4"
    outer = os.environ.get("CUMCM_OUTERN", None)

    r = one(port, seed, nsrc, direc, mode, outer=outer)
    if not r:
        return
    real = r["real"]
    real_ch = set(real)
    missed = real_ch - r["cleared"]
    print("=== %s seed=%d sources=%d directional=%.2f outer_n=%s ===" % (mode, seed, nsrc, direc, os.environ.get("CUMCM_OUTERN", "6")))
    print("虚拟时间         : %s" % r["vt"])
    print("真实源数         : %d   (%d 个定向)" %
          (len(real), sum(1 for s in real.values() if s["directional"])))
    print("粗扫检测到信号   : %d  %s" % (len(r["found"]), sorted(r["found"])))
    print("清除成功         : %d  %s" % (len(r["cleared"]), sorted(r["cleared"])))
    print("★ 漏源           : %d  %s" % (len(missed), sorted(missed)))
    print()
    if missed:
        print("漏源明细：")
        for c in sorted(missed):
            s = real[c]
            got = "收到过信号但没清" if c in r["found"] else "从未收到信号"
            print("  ch%-3d 定向=%-6s recv_r=%6.0f  r=%6.0f  位置(%7.0f,%7.0f)  %s"
                  % (c, s["directional"], s["recv_r"],
                     math.hypot(s["pos"][0], s["pos"][1]),
                     s["pos"][0], s["pos"][1], got))
    # 误报：清了不存在的频道（不该发生，但值得查）
    ghost = r["cleared"] - real_ch
    if ghost:
        print("\n!! 清除了不存在的频道:", sorted(ghost))


if __name__ == "__main__":
    main()
