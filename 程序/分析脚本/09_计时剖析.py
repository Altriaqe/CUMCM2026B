"""时间分解：把虚拟时间按「扫描移动 / 扫描测量 / 清除移动 / 清除测量 / 清除动作」归因。"""
import threading, math, sys
import mock_simulator as ms
from client import RobotClient
import strategy as st

class PClient(RobotClient):
    """包一层，记录每次调用引起的虚拟时间增量。"""
    def __init__(self, rid):
        super().__init__(rid)
        self.buckets = {}
        self.cur = "init"
        self._last_vt = 0.0
        self.pos = (0.0, 0.0)

    def _acc(self, tag, r):
        vt = r.get("virtual_time_s", self._last_vt)
        d = vt - self._last_vt
        if d > 0:
            self.buckets[tag] = self.buckets.get(tag, 0.0) + d
        self._last_vt = vt
        return r

    def measure(self, x, y, ch):
        tag = "scan_move" if math.dist(self.pos, (x, y)) > 1.0 else "scan_stay"
        r = super().measure(x, y, ch)
        self.pos = (x, y)
        return self._acc("scan:" + tag + ":meas", r)

    def clear(self, x, y, ch):
        r = super().clear(x, y, ch)
        self.pos = (x, y)
        return self._acc("clear_act", r)


def serve(port, n, dr, seed):
    import random
    random.seed(seed)
    ms.Handler.world = ms.World(n, dr, seed)
    ms.Handler.lock = threading.Lock()
    srv = ms.ThreadingHTTPServer(("127.0.0.1", port), ms.Handler)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    return srv


def run(seed, n, dr, port):
    srv = serve(port, n, dr, seed)
    c = PClient("test")
    c.enter()
    pts = st.coverage_points()
    sightings, _ = st.two_pass_sweep(c, pts)
    est = st.estimate_sources(sightings)
    order = st.clear_order((0.0, 0.0), est)
    cleared = 0
    for ch, x, y in order:
        if st.homing_clear(c, ch, (x, y)):
            cleared += 1
    vt = c.virtual_time
    c.exit()
    srv.shutdown()
    return vt, cleared, c.buckets, len(est)


if __name__ == "__main__":
    seeds = [1, 2, 3, 4, 5, 6, 7, 8]
    port = 8850
    tot = {}
    vts, cls = [], []
    for sd in seeds:
        port += 1
        vt, cl, b, ne = run(sd, 13, 0.4, port)
        vts.append(vt); cls.append(cl)
        for k, v in b.items():
            tot[k] = tot.get(k, 0.0) + v
        print(f"seed{sd}: vt={vt:.0f} cleared={cl} est={ne}")
    print()
    print("=== 汇总（%d 个种子，共 %d 个源）===" % (len(seeds), sum(cls)))
    s = sum(tot.values())
    for k in sorted(tot, key=lambda k: -tot[k]):
        print(f"  {k:24s} {tot[k]:9.1f} s   {100*tot[k]/s:5.1f}%")
    print(f"  {'TOTAL':24s} {s:9.1f} s")
    print(f"  平均每源 {s/sum(cls):.0f} s")
