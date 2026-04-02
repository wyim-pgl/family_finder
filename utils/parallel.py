"""Multiprocessing wrapper for parallel orthogroup processing."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, Dict, List, Tuple
import logging

logger = logging.getLogger("family_finder")


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

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(func, item): i for i, item in enumerate(items)}
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                idx = futures[future]
                logger.error(f"Worker failed on item {idx}: {e}")
                results.append(None)

    return results
