from __future__ import annotations

from pathlib import Path


SCRIPT = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/patel_2021_3hfm/run_patellike_3hfm_watcher.sh"
)


def test_watcher_wrapper_uses_conservative_defaults() -> None:
    text = SCRIPT.read_text(encoding="utf-8")

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-13}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-2}"' in text
    assert 'MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"' in text
    assert 'SKIP_CHARGE_CHANGING="${SKIP_CHARGE_CHANGING:-1}"' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert '--max-launches-per-pass "${MAX_LAUNCHES_PER_PASS}"' in text
    assert '--skip-charge-changing' in text
    assert 'watch_patellike_3hfm.py' in text
