"""
本地 mock 模拟器：实现 B 题 HTTP 协议，供赛前联调策略，不依赖官方模拟器。

用法：
    python mock_simulator.py [--port 2026] [--seed 1] [--sources 13] [--directional 0.0]

然后另开终端运行 robot.py（ROBOT_ID 填任意非空串，mock 不校验队号）。
mock 会打印真实干扰源位置，供复盘对照；正式比赛这些真值不可见。
"""
import argparse
import hashlib
import json
import math
import random
import struct
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

AREA_R = 1800.0
SPEED = 5.0
NEAR_R = 5.0        # 距离 <= 5m 且在覆盖锥内 -> near
CLEAR_R = 20.0      # 清除半径


def _det_unit(*vals):
    """确定性伪随机数：同样的输入永远得到同样的 [0,1) 输出。

    用于模拟题设中「同一地点重测误差不变」的示向度误差 —— 误差是位置的
    函数，不是每次调用的随机噪声。"""
    h = hashlib.md5()
    for v in vals:
        h.update(struct.pack("<d", float(v)))   # 定点化，避免浮点尾数抖动
    return int.from_bytes(h.digest()[:8], "big") / 2.0 ** 64


class Source:
    def __init__(self, channel, pos, directional, seed):
        self.channel = channel
        self.pos = pos
        self.directional = directional
        self.recv_r = random.uniform(1000, 1500)
        self.direction = random.uniform(0, 360) if directional else None
        self.cleared = False
        # 每个源一个固定偏移：误差随位置连续变化，但同一位置必然相同
        self.err_phase = _det_unit(seed, channel, 1)
        # 幅值在 [0.6, 1.0] 宽带内取值，保证 ±1° 是实的上界而非摆设
        self.err_gain = 0.6 + 0.4 * _det_unit(seed, channel, 2)

    # 误差空间尺度：500m。误差在 ~500m 尺度上缓慢起伏，几米内几乎不变，
    # 既满足「同地点重测不变」，又使相距几百米的示向度携带独立信息。
    ERR_SCALE = 500.0

    def err_at(self, x, y):
        """检测点 (x,y) 处的示向度误差，单位度，|err| <= 1。

        误差是 (检测点, 源) 的确定性连续函数（值噪声 + 平滑插值）：
        同地点重测结果逐位相同，且相邻位置之间连续无跳变。"""
        g = self.ERR_SCALE
        fx, fy = x / g, y / g
        ix, iy = math.floor(fx), math.floor(fy)
        tx, ty = fx - ix, fy - iy
        sx = tx * tx * (3.0 - 2.0 * tx)      # smoothstep
        sy = ty * ty * (3.0 - 2.0 * ty)
        c = self.channel * 1.0
        h = [_det_unit(c, self.err_phase, ix + dx, iy + dy) * 2.0 - 1.0
             for dx, dy in ((0, 0), (1, 0), (0, 1), (1, 1))]
        v = ((1 - sx) * (1 - sy) * h[0] + sx * (1 - sy) * h[1] +
             (1 - sx) * sy * h[2] + sx * sy * h[3])
        return self.err_gain * v


class World:
    def __init__(self, n_sources, directional_ratio, seed):
        channels = random.sample(range(1, 21), n_sources)
        self.sources = []
        for c in channels:
            self.sources.append(Source(c, _rand_in_circle(AREA_R),
                                       random.random() < directional_ratio,
                                       seed))
        self.robot = (0.0, 0.0)
        self.channel = 1
        self.vt = 0.0
        self.entered = False


def _rand_in_circle(r):
    while True:
        x = random.uniform(-r, r)
        y = random.uniform(-r, r)
        if x * x + y * y <= r * r:
            return (x, y)


def _source_on(world, channel):
    for s in world.sources:
        if s.channel == channel and not s.cleared:
            return s
    return None


def _in_cone(s, pos):
    if not s.directional:
        return True
    ang = math.degrees(math.atan2(pos[1] - s.pos[1], pos[0] - s.pos[0])) % 360
    diff = (ang - s.direction) % 360
    return diff <= 90 or diff >= 270


def _move_cost(world, pos):
    return math.hypot(pos[0] - world.robot[0], pos[1] - world.robot[1]) / SPEED


def _reject():
    return {"accepted": False, "real_timestamp_ms": 0, "virtual_time_s": 0}


def _enter(world):
    world.entered = True
    return {"accepted": True, "real_timestamp_ms": 0, "virtual_time_s": 0.0,
            "max_virtual_duration_s": 360000, "max_real_duration_s": 1200,
            "remaining_real_duration_s": 1200}


def _measure(world, x, y, channel):
    if not world.entered:
        return _reject()
    world.vt += _move_cost(world, (x, y)) + (1.0 if channel != world.channel else 0.0) + 5.0
    world.robot = (x, y)
    world.channel = channel
    r = {"accepted": True, "real_timestamp_ms": 0,
         "virtual_time_s": round(world.vt, 6)}
    s = _source_on(world, channel)
    if s is None:
        r["measure_result"] = "no_signal"
    else:
        d = math.hypot(x - s.pos[0], y - s.pos[1])
        if d <= NEAR_R and _in_cone(s, (x, y)):
            r["measure_result"] = "near"
        elif d <= s.recv_r and _in_cone(s, (x, y)):
            true = math.degrees(math.atan2(s.pos[1] - y, s.pos[0] - x)) % 360
            r["measure_result"] = "direction"
            r["svd_deg"] = round((true + s.err_at(x, y)) % 360, 2)
        else:
            r["measure_result"] = "no_signal"
    return r


def _clear(world, x, y, channel):
    if not world.entered:
        return _reject()
    world.vt += _move_cost(world, (x, y))
    s = _source_on(world, channel)
    ok = s is not None and math.hypot(x - s.pos[0], y - s.pos[1]) <= CLEAR_R
    if ok:
        s.cleared = True
    world.vt += 5.0 if ok else 3.0
    world.robot = (x, y)
    return {"accepted": True, "real_timestamp_ms": 0,
            "virtual_time_s": round(world.vt, 6),
            "clear_result": "success" if ok else "no_target_in_range"}


def _exit(world):
    return {"accepted": True, "real_timestamp_ms": 0, "virtual_time_s": world.vt,
            "exit_reason": "user_exit"}


class Handler(BaseHTTPRequestHandler):
    world = None

    def _send(self, obj, code=200):
        body = json.dumps(obj).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
        except Exception:
            return self._send(_reject(), 400)
        w = Handler.world
        path = self.path
        # 世界状态可变，且 mock 用 ThreadingHTTPServer，必须串行化
        with Handler.lock:
            if path == "/enter":
                return self._send(_enter(w))
            p = payload.get("position", {})
            x = p.get("x")
            y = p.get("y")
            c = payload.get("channel")
            if path == "/measure":
                if x is None or y is None or c is None:
                    return self._send(_reject(), 400)
                return self._send(_measure(w, x, y, c))
            if path == "/clear":
                if x is None or y is None or c is None:
                    return self._send(_reject(), 400)
                return self._send(_clear(w, x, y, c))
            if path == "/exit":
                return self._send(_exit(w))
            return self._send(_reject(), 404)

    def log_message(self, *a):
        pass

    def do_GET(self):
        # 便于本地查看世界状态（官方模拟器无此接口）
        if self.path == "/_truth":
            w = Handler.world
            return self._send({"sources": [
                {"channel": s.channel, "pos": s.pos, "directional": s.directional,
                 "direction": s.direction, "recv_r": s.recv_r,
                 "cleared": s.cleared} for s in w.sources],
                "virtual_time_s": w.vt, "robot": w.robot})
        return self._send(_reject(), 404)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=2026)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sources", type=int, default=13)
    ap.add_argument("--directional", type=float, default=0.0)
    args = ap.parse_args()
    random.seed(args.seed)
    Handler.world = World(args.sources, args.directional, args.seed)
    Handler.lock = threading.Lock()
    print(f"mock 模拟器启动，端口 {args.port}，干扰源 {args.sources} 个，"
          f"定向比例 {args.directional}，种子 {args.seed}")
    print("真实干扰源位置（仅供本地复盘对照）：")
    for s in Handler.world.sources:
        extra = (f"定向 dir={s.direction:.0f}" if s.directional else "全向")
        print(f"  频道{s.channel}: ({s.pos[0]:.1f}, {s.pos[1]:.1f}) "
              f"recv={s.recv_r:.0f} {extra}")
    ThreadingHTTPServer(("127.0.0.1", args.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
