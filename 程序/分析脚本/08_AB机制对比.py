"""A/B：MEC 停机判据 vs 原版 homing_clear。"""
import threading, sys, importlib, io, contextlib
import mock_simulator as ms
from client import RobotClient
import strategy as st

LEGACY = st.homing_clear
def _legacy(*a, **kw):
    kw.pop("stop_radius", None)
    return LEGACY(*a, **kw)

def serve(port, n, dr, seed):
    import random
    random.seed(seed)
    ms.Handler.world = ms.World(n, dr, seed)
    ms.Handler.lock = threading.Lock()
    srv = ms.ThreadingHTTPServer(("127.0.0.1", port), ms.Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv

def run(seed, fn, port, n=13, dr=0.4):
    st.homing_clear = fn
    srv = serve(port, n, dr, seed)
    c = RobotClient("test"); c.enter()
    pts = st.coverage_points()
    sig, _ = st.two_pass_sweep(c, pts)
    est = st.estimate_sources(sig)
    order = st.clear_order((0.0, 0.0), est)
    cl = 0
    for ch, x, y in order:
        if st.homing_clear(c, ch, (x, y)):
            cl += 1
    vt = c.virtual_time; c.exit(); srv.shutdown()
    return vt, cl, len(est)

if __name__ == "__main__":
    seeds = [1, 2, 3, 4, 5, 6]
    mode = sys.argv[1]
    fn = st.homing_clear if mode == "new" else _legacy
    tot_v = tot_c = 0
    port = 8900 + int(sys.argv[2]) * 50
    for sd in seeds:
        port += 1
        vt, cl, ne = run(sd, fn, port)
        tot_v += vt; tot_c += cl
        print("  seed%d: vt=%7.1f cleared=%2d/%d" % (sd, vt, cl, ne), flush=True)
    print("%s: 总虚拟时间 %.0f s, 清除 %d, 平均 %.1f s/源"
          % (mode, tot_v, tot_c, tot_v / max(tot_c, 1)))
