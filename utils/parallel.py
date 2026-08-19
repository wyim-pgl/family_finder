"""Multiprocessing wrapper for parallel orthogroup processing.

Worker logging (issue #19): logging handlers configured in the parent are NOT
inherited by ProcessPoolExecutor children, so every per-OG diagnostic
(pruning results, stop-codon removals, TreeShrink warnings) used to be
silently discarded — pipeline.log from the 5sp run contains zero worker
lines despite --verbose. parallel_map now wires a QueueHandler in each
worker to a QueueListener in the parent, so worker log records reach the
real handlers.
"""

import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from logging.handlers import QueueHandler, QueueListener
from typing import Any, Callable, List, Tuple

logger = logging.getLogger("family_finder")


def _worker_log_init(queue, level: int) -> None:
    """Executor initializer: route the worker's family_finder logger into the
    parent's queue. Replaces any handlers the child inherited (fork) or lacks
    (spawn) so records are neither duplicated nor lost."""
    wlogger = logging.getLogger("family_finder")
    wlogger.handlers = [QueueHandler(queue)]
    wlogger.setLevel(level)
    wlogger.propagate = False


def parallel_map(
    func: Callable,
    items: List[Tuple],
    n_workers: int = 8,
) -> List[Any]:
    """Map a function over items in parallel using ProcessPoolExecutor.

    Args:
        func: Function to call. Must accept a single tuple of args (picklable).
        items: List of argument tuples to pass to func.
        n_workers: Number of parallel workers.

    Returns:
        List of results in completion order.
    """
    results = []

    if n_workers <= 1 or len(items) <= 1:
        for item in items:
            results.append(func(item))
        return results

    # Relay worker log records to the parent's handlers (fall back to the
    # root logger's handlers when family_finder has none of its own).
    handlers = logger.handlers or logging.getLogger().handlers
    log_queue = multiprocessing.Manager().Queue(-1)
    listener = QueueListener(log_queue, *handlers, respect_handler_level=True)
    listener.start()

    try:
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_worker_log_init,
            initargs=(log_queue, logger.getEffectiveLevel()),
        ) as executor:
            futures = {executor.submit(func, item): i for i, item in enumerate(items)}
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as e:
                    idx = futures[future]
                    logger.error(f"Worker failed on item {idx}: {e}")
                    results.append(None)
    finally:
        listener.stop()

    return results
