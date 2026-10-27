# -*- coding: utf-8 -*-
"""A/B：同一 mock 世界（同 seed）下，顺路分流开/关的总虚拟时间与路程对比。

用法：
    python _att/ab_divert.py <port> <seed> [sources] [directional]

前提：mock 模拟器已在该 port 上以该 seed 启动（本脚本会提示重启）。
因为 mock 的 World 在进程启动时构造，脚本无法重置，所以每一轮都必须
重启 mock。这里直接用 subprocess 起停，保证干净。
"""
import json, os, re, socket, subprocess, sys, time
sys.stdout = __import__('io').TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable


def wait_port(port, timeout=15.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        s = socket.socket()
        s.settimeout(0.4)
        try:
            s.connect(("127.0.0.1", port))
            s.close()
            return True
        except OSError:
            time.sleep(0.3)
        finally:
            try:
                s.close()
            except Exception:
                pass
    return False


def kill_port(port):
    try:
        out = subprocess.run(["netstat", "-ano"], capture_output=True, text=True,
                             timeout=20).stdout
    except Exception:
        return
    pids = set()
    for line in out.splitlines():
        if "127.0.0.1:%d " % port in line and "LISTENING" in line:
            pids.add(line.split()[-1])
    for p in pids:
        subprocess.run(["taskkill", "/F", "/PID", p], capture_output=True)


def run_once(port, seed, sources, directional, divert, tag):
    kill_port(port)
    time.sleep(0.5)
    env = dict(os.environ)
    env["CUMCM_DIVERT"] = "1" if divert else "0"
    mock = subprocess.Popen([PY, "-u", "mock_simulator.py", "--port", str(port),
                             "--seed", str(seed), "--sources", str(sources),
                             "--directional", str(directional)],
                            cwd=ROOT, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    if not wait_port(port):
        print("mock 起不来")
        return None
    try:
        t0 = time.time()
        pr = subprocess.run([PY, "robot.py", "--id", tag, "--port", str(port)],
                            cwd=ROOT, capture_output=True, text=True,
                            timeout=1800, env=env)
        text = pr.stdout
    finally:
        mock.terminate()
        try:
            mock.wait(timeout=5)
        except Exception:
            mock.kill()
    lines = {}
    for ln in text.splitlines():
        for key in ("退出", "总路程", "平均定位"):
            if key in ln:
                lines[key] = ln.strip()
    m = re.search(r"虚拟时间: ([\d.]+) \| 清除: (\d+)", text)
    dist = re.search(r"总路程: (\d+)m", text)
    mov = re.search(r"移动: ([\d.]+)s", text)
    meas = re.search(r"测量×(\d+)", text)
    return {
        "vt": float(m.group(1)) if m else None,
        "cleared": int(m.group(2)) if m else None,
        "dist": int(dist.group(1)) if dist else None,
        "move": float(mov.group(1)) if mov else None,
        "meas": int(meas.group(1)) if meas else None,
        "wall": time.time() - t0,
    }


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 2041
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 42
    nsrc = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    direc = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0

    print("seed=%d sources=%d directional=%.2f" % (seed, nsrc, direc))
    print("%-10s %10s %8s %8s %8s %8s %8s" %
          ("模式", "虚拟时间", "清除", "路程", "移动", "测量", "用时"))
    res = {}
    for divert, tag in ((False, "OFF"), (True, "ON")):
        r = run_once(port, seed, nsrc, direc, divert, "ab" + tag)
        res[tag] = r
        if r:
            print("%-10s %10.1f %8s %8s %8.1f %8s %8.0f" %
                  (tag, r["vt"], r["cleared"], r["dist"], r["move"], r["meas"],
                   r["wall"]))
    if res.get("OFF") and res.get("ON"):
        a, b = res["OFF"], res["ON"]
        print()
        print("分流 ON 相对 OFF:  虚拟时间 %+.1f (%.1f%%)   总路程 %+d (%.1f%%)"
              % (b["vt"] - a["vt"], 100 * (b["vt"] - a["vt"]) / a["vt"],
                 b["dist"] - a["dist"], 100 * (b["dist"] - a["dist"]) / a["dist"]))


if __name__ == "__main__":
    main()
