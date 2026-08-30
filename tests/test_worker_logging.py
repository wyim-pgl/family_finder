"""Issue #19: worker log records must reach the parent's handlers."""

import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import utils.parallel as parallel  # noqa: E402
from utils.parallel import parallel_map  # noqa: E402


def _noisy_worker(item):
    # Arrange (in worker): module-level so it is picklable
    logging.getLogger("family_finder").warning(f"worker-log item={item[0]}")
    return item[0] * 2


class _Capture(logging.Handler):
    def __init__(self):
        super().__init__()
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


def test_worker_logs_reach_parent_handlers():
    # Arrange
    lg = logging.getLogger("family_finder")
    cap = _Capture()
    old_handlers, old_level = lg.handlers[:], lg.level
    lg.handlers = [cap]
    lg.setLevel(logging.INFO)
    try:
        # Act: real ProcessPoolExecutor round-trip (2 workers, 3 items)
        results = parallel_map(_noisy_worker, [(1,), (2,), (3,)], n_workers=2)
    finally:
        lg.handlers, lg.level = old_handlers, old_level

    # Assert
    assert sorted(results) == [2, 4, 6]
    got = sorted(m for m in cap.messages if m.startswith("worker-log"))
    assert got == ["worker-log item=1", "worker-log item=2", "worker-log item=3"]


def test_serial_path_still_logs_directly():
    lg = logging.getLogger("family_finder")
    cap = _Capture()
    old_handlers, old_level = lg.handlers[:], lg.level
    lg.handlers = [cap]
    lg.setLevel(logging.INFO)
    try:
        results = parallel_map(_noisy_worker, [(7,)], n_workers=8)
    finally:
        lg.handlers, lg.level = old_handlers, old_level
    assert results == [14]
    assert "worker-log item=7" in cap.messages


def test_manager_eof_falls_back_to_plain_queue_and_worker_logs_still_reach_parent(
        monkeypatch):
    # Arrange: in a restricted sandbox the Manager child gets PermissionError
    # while opening its listener socket, which reaches the parent as EOFError.
    # The fallback must handle the parent-side exception, not only the child's
    # underlying error type.
    def manager_socket_was_refused():
        raise EOFError("manager child exited before sending its address")

    monkeypatch.setattr(parallel.multiprocessing, "Manager",
                        manager_socket_was_refused)
    lg = logging.getLogger("family_finder")
    cap = _Capture()
    old_handlers, old_level = lg.handlers[:], lg.level
    lg.handlers = [cap]
    lg.setLevel(logging.INFO)
    try:
        # Act: use the real plain multiprocessing.Queue and executor after the
        # forced Manager failure, so this proves the complete fallback path.
        results = parallel_map(_noisy_worker, [(1,), (2,)], n_workers=2)
    finally:
        lg.handlers, lg.level = old_handlers, old_level

    # Assert
    assert sorted(results) == [2, 4]
    got = sorted(m for m in cap.messages if m.startswith("worker-log"))
    assert got == ["worker-log item=1", "worker-log item=2"]
