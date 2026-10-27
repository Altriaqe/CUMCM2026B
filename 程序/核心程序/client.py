"""
模拟器 HTTP 客户端。

封装 /enter /measure /clear /exit 四条指令，处理幂等重试、日志与虚拟时间记录。

协议要点（详见 B 题附件 2）：
- 每个新动作使用新的 request_id；仅在网络异常重试同一动作时复用原 request_id 与原请求体。
- 必须同时检查 HTTP 状态码与响应体 accepted 字段。
- accepted=false 时 virtual_time_s 为 0，不代表当前虚拟时刻。
- 只有成功的 /measure 会更新测向机当前频道；/clear 不更新。
"""
import json
import time
import urllib.request
import urllib.error


class NoResponse(Exception):
    """网络层失败：没拿到 HTTP 响应（连接被关、超时、接口未开放等）。"""


class SimError(Exception):
    """业务层失败：拿到了响应，但 HTTP 非 200 或 accepted=false。"""


class RobotClient:
    def __init__(self, robot_id, base_url="http://127.0.0.1:2026", timeout=5.0,
                 trace_path=None):
        self.robot_id = robot_id
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self._seq = 0
        self.log = []          # [(path, payload, response_dict), ...]
        self.virtual_time = 0.0
        self.position = (0.0, 0.0)
        self.current_channel = 1
        # 统计计数
        self.n_measure = 0
        self.n_clear_ok = 0
        self.n_clear_fail = 0
        self.n_channel_switch = 0
        self.total_dist = 0.0
        self._last_pos = (0.0, 0.0)
        # 逐动作轨迹落盘：官方加密日志无法解读，只能自己记。
        # 每次动作后立即 flush，程序被强杀也不丢已经发生的动作。
        self.trace_path = trace_path
        self.trace = None
        self.phase = ""        # 当前阶段标签，由调用方通过 note() 设置
        if trace_path:
            self.trace = open(trace_path, "w", encoding="utf-8")

    def note(self, tag):
        """设置阶段标签，会写进此后每条轨迹记录，便于事后定位 fail 发生在哪一段。"""
        self.phase = str(tag)
        if self.trace:
            self.trace.write(json.dumps({"kind": "mark", "tag": self.phase},
                                        ensure_ascii=False) + "\n")
            self.trace.flush()

    def _trace(self, path, payload, resp):
        if not self.trace:
            return
        p = payload.get("position") or {}
        rec = {
            "kind": "act",
            "path": path,
            "rid": payload.get("request_id"),
            "vt": resp.get("virtual_time_s"),
            "phase": self.phase,
            "x": p.get("x"), "y": p.get("y"),
            "ch": payload.get("channel"),
            "res": resp.get("measure_result") or resp.get("clear_result")
                   or resp.get("exit_reason"),
            "svd": resp.get("svd_deg"),
            "acc": resp.get("accepted"),
            "pos": list(self.position),
        }
        try:
            self.trace.write(json.dumps(rec, ensure_ascii=False) + "\n")
            self.trace.flush()
        except Exception:
            pass               # 轨迹写失败绝不能影响比赛流程

    def close(self):
        if self.trace:
            try:
                self.trace.close()
            except Exception:
                pass
            self.trace = None

    # ---- 底层 ----
    def _next_id(self, tag):
        self._seq += 1
        return f"{tag}-{self._seq}"

    def _base(self, request_id):
        return {"arena_id": "default", "robot_id": self.robot_id,
                "request_id": request_id}

    def _post(self, path, payload):
        data = json.dumps(payload, separators=(",", ":")).encode("utf-8")
        req = urllib.request.Request(
            self.base_url + path, data=data,
            headers={"Content-Type": "application/json"}, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                return r.status, json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            try:
                return e.code, json.loads(e.read().decode("utf-8"))
            except Exception:
                return e.code, {}
        except Exception as e:          # URLError / 超时 / 连接被直接关闭
            raise NoResponse(f"{path}: {e}") from e

    def _call(self, path, payload, request_id):
        """发送动作；仅在网络层失败时以相同 request_id 重试（保证幂等）。"""
        while True:
            try:
                status, resp = self._post(path, payload)
            except NoResponse:
                time.sleep(0.5)
                continue
            self.log.append((path, payload, resp))
            self._trace(path, payload, resp)
            return status, resp

    def _ok(self, status, resp):
        if status != 200:
            raise SimError(f"HTTP {status}: {resp}")
        if not resp.get("accepted"):
            raise SimError(f"accepted=false: {resp}")
        if resp.get("virtual_time_s"):
            self.virtual_time = resp["virtual_time_s"]
        return resp

    # ---- 动作 ----
    def enter(self):
        rid = self._next_id("enter")
        r = self._ok(*self._call("/enter", self._base(rid), rid))
        self.position = (0.0, 0.0)
        self.current_channel = 1
        return r

    def measure(self, x, y, channel):
        rid = self._next_id("measure")
        p = self._base(rid)
        p["position"] = {"x": x, "y": y}
        p["channel"] = channel
        r = self._ok(*self._call("/measure", p, rid))
        import math
        self.total_dist += math.hypot(x - self._last_pos[0], y - self._last_pos[1])
        if channel != self.current_channel:
            self.n_channel_switch += 1
        self.n_measure += 1
        self.position = (x, y)
        self._last_pos = (x, y)
        self.current_channel = channel
        return r

    def clear(self, x, y, channel):
        rid = self._next_id("clear")
        p = self._base(rid)
        p["position"] = {"x": x, "y": y}
        p["channel"] = channel
        r = self._ok(*self._call("/clear", p, rid))
        import math
        self.total_dist += math.hypot(x - self._last_pos[0], y - self._last_pos[1])
        if r.get("clear_result") == "success":
            self.n_clear_ok += 1
        else:
            self.n_clear_fail += 1
        self.position = (x, y)
        self._last_pos = (x, y)
        return r

    def exit(self):
        rid = self._next_id("exit")
        return self._ok(*self._call("/exit", self._base(rid), rid))
