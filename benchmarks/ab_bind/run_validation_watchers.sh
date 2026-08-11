#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
SPLIT_FILE="${ROOT}/benchmarks/ab_bind/splits/ab_bind_rbfe_core_v1_split_v1.yml"
ROBUST_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"
PRIORITY_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"
TARGETED_REPEAT_SPREAD_RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"
TARGETED_LAMBDA_RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"
SAMPLING_QC_RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"
DEEP_RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"
ULTRA_RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"
STALE_RECOVERY_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
GAP_RECOVERY_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
BACKLOG_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"

ROBUST_LOG="${ROBUST_ROOT}/reports/watch/watch_validation_priority_merged_auto.log"
PRIORITY_LOG="${PRIORITY_ROOT}/reports/watch/watch_validation_priority_relaunched_auto.log"
RESCUE_LOG="${RESCUE_ROOT}/reports/watch/watch_validation_priority_rescue_auto.log"
TARGETED_REPEAT_SPREAD_RESCUE_LOG="${TARGETED_REPEAT_SPREAD_RESCUE_ROOT}/reports/watch/watch_validation_targeted_repeat_spread_rescue_auto.log"
TARGETED_LAMBDA_RESCUE_LOG="${TARGETED_LAMBDA_RESCUE_ROOT}/reports/watch/watch_validation_targeted_lambda_rescue_auto.log"
SAMPLING_QC_RESCUE_LOG="${SAMPLING_QC_RESCUE_ROOT}/reports/watch/watch_validation_sampling_qc_rescue_auto.log"
DEEP_RESCUE_LOG="${DEEP_RESCUE_ROOT}/reports/watch/watch_validation_deep_rescue_auto.log"
ULTRA_RESCUE_LOG="${ULTRA_RESCUE_ROOT}/reports/watch/watch_validation_ultra_rescue_auto.log"
STALE_RECOVERY_LOG="${STALE_RECOVERY_ROOT}/reports/watch/watch_validation_stale_recovery_auto.log"
GAP_RECOVERY_LOG="${GAP_RECOVERY_ROOT}/reports/watch/watch_validation_gap_recovery_auto.log"
BACKLOG_LOG="${BACKLOG_ROOT}/reports/watch/watch_validation_priority_backlog_auto.log"

ROBUST_PIDFILE="${ROBUST_ROOT}/reports/watch/watch_validation_robust.pid"
PRIORITY_PIDFILE="${PRIORITY_ROOT}/reports/watch/watch_validation_priority.pid"
RESCUE_PIDFILE="${RESCUE_ROOT}/reports/watch/watch_validation_rescue.pid"
TARGETED_REPEAT_SPREAD_RESCUE_PIDFILE="${TARGETED_REPEAT_SPREAD_RESCUE_ROOT}/reports/watch/watch_validation_targeted_repeat_spread_rescue.pid"
TARGETED_LAMBDA_RESCUE_PIDFILE="${TARGETED_LAMBDA_RESCUE_ROOT}/reports/watch/watch_validation_targeted_lambda_rescue.pid"
SAMPLING_QC_RESCUE_PIDFILE="${SAMPLING_QC_RESCUE_ROOT}/reports/watch/watch_validation_sampling_qc_rescue.pid"
DEEP_RESCUE_PIDFILE="${DEEP_RESCUE_ROOT}/reports/watch/watch_validation_deep_rescue.pid"
ULTRA_RESCUE_PIDFILE="${ULTRA_RESCUE_ROOT}/reports/watch/watch_validation_ultra_rescue.pid"
STALE_RECOVERY_PIDFILE="${STALE_RECOVERY_ROOT}/reports/watch/watch_validation_stale_recovery.pid"
GAP_RECOVERY_PIDFILE="${GAP_RECOVERY_ROOT}/reports/watch/watch_validation_gap_recovery.pid"
BACKLOG_PIDFILE="${BACKLOG_ROOT}/reports/watch/watch_validation_priority_backlog.pid"

ROBUST_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_robust_watcher.sh"
PRIORITY_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_priority_watcher.sh"
RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_rescue_watcher.sh"
TARGETED_REPEAT_SPREAD_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_targeted_repeat_spread_rescue_watcher.sh"
TARGETED_LAMBDA_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_targeted_lambda_rescue_watcher.sh"
SAMPLING_QC_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_sampling_qc_rescue_watcher.sh"
DEEP_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_deep_rescue_watcher.sh"
ULTRA_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_ultra_rescue_watcher.sh"
STALE_RECOVERY_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_stale_recovery_watcher.sh"
GAP_RECOVERY_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_gap_recovery_watcher.sh"
BACKLOG_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_priority_backlog_watcher.sh"

ROBUST_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${ROBUST_ROOT}"
PRIORITY_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${PRIORITY_ROOT}"
RESCUE_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${RESCUE_ROOT}"
TARGETED_REPEAT_SPREAD_RESCUE_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${TARGETED_REPEAT_SPREAD_RESCUE_ROOT}"
TARGETED_LAMBDA_RESCUE_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${TARGETED_LAMBDA_RESCUE_ROOT}"
SAMPLING_QC_RESCUE_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${SAMPLING_QC_RESCUE_ROOT}"
DEEP_RESCUE_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${DEEP_RESCUE_ROOT}"
ULTRA_RESCUE_WATCHER_PATTERN="watch_validation_priority.py --plan-root ${ULTRA_RESCUE_ROOT}"
STALE_RECOVERY_WATCHER_PATTERN="watch_validation_priority.py .*--split-file ${SPLIT_FILE} .*--watch-tag stale"
GAP_RECOVERY_WATCHER_PATTERN="watch_validation_priority.py .*--split-file ${SPLIT_FILE} .*--watch-tag gap"
BACKLOG_WATCHER_PATTERN="watch_validation_priority.py .*--split-file ${SPLIT_FILE} .*--watch-tag backlog"

ACTION="${1:-status}"
if [ "$#" -gt 0 ]; then
  shift
fi

ensure_watch_dirs() {
  mkdir -p \
    "$(dirname "${ROBUST_LOG}")" \
    "$(dirname "${PRIORITY_LOG}")" \
    "$(dirname "${RESCUE_LOG}")" \
    "$(dirname "${TARGETED_REPEAT_SPREAD_RESCUE_LOG}")" \
    "$(dirname "${TARGETED_LAMBDA_RESCUE_LOG}")" \
    "$(dirname "${SAMPLING_QC_RESCUE_LOG}")" \
    "$(dirname "${DEEP_RESCUE_LOG}")" \
    "$(dirname "${ULTRA_RESCUE_LOG}")" \
    "$(dirname "${STALE_RECOVERY_LOG}")" \
    "$(dirname "${GAP_RECOVERY_LOG}")" \
    "$(dirname "${BACKLOG_LOG}")"
}

pid_running() {
  local pid="$1"
  [ -n "${pid}" ] && kill -0 "${pid}" 2>/dev/null
}

watcher_pattern_for_target() {
  local target="$1"
  case "${target}" in
    robust)
      printf '%s\n' "${ROBUST_WATCHER_PATTERN}"
      ;;
    priority)
      printf '%s\n' "${PRIORITY_WATCHER_PATTERN}"
      ;;
    rescue)
      printf '%s\n' "${RESCUE_WATCHER_PATTERN}"
      ;;
    targeted)
      printf '%s\n' "${TARGETED_REPEAT_SPREAD_RESCUE_WATCHER_PATTERN}"
      ;;
    lambda)
      printf '%s\n' "${TARGETED_LAMBDA_RESCUE_WATCHER_PATTERN}"
      ;;
    sampling-qc)
      printf '%s\n' "${SAMPLING_QC_RESCUE_WATCHER_PATTERN}"
      ;;
    deep)
      printf '%s\n' "${DEEP_RESCUE_WATCHER_PATTERN}"
      ;;
    ultra)
      printf '%s\n' "${ULTRA_RESCUE_WATCHER_PATTERN}"
      ;;
    stale)
      printf '%s\n' "${STALE_RECOVERY_WATCHER_PATTERN}"
      ;;
    gap)
      printf '%s\n' "${GAP_RECOVERY_WATCHER_PATTERN}"
      ;;
    backlog)
      printf '%s\n' "${BACKLOG_WATCHER_PATTERN}"
      ;;
  esac
}

find_running_pid_by_patterns() {
  local pattern
  local pid
  for pattern in "$@"; do
    if [ -z "${pattern}" ]; then
      continue
    fi
    pid="$(pgrep -f "${pattern}" | head -n 1 || true)"
    if pid_running "${pid}"; then
      printf '%s\n' "${pid}"
      return 0
    fi
  done
  return 1
}

status_line() {
  local label="$1"
  local pidfile="$2"
  local pattern="$3"
  local watcher_pattern
  watcher_pattern="$(watcher_pattern_for_target "${label}")"
  if [ -f "${pidfile}" ]; then
    local pid
    pid="$(cat "${pidfile}")"
    if pid_running "${pid}"; then
      printf '%s\tpid=%s\tstatus=running\n' "${label}" "${pid}"
      return
    fi
  fi
  local pid
  pid="$(find_running_pid_by_patterns "${pattern}" "${watcher_pattern}" || true)"
  if pid_running "${pid}"; then
    printf '%s\tpid=%s\tstatus=running(no-pidfile)\n' "${label}" "${pid}"
  else
    printf '%s\tstatus=stopped\n' "${label}"
  fi
}

start_one() {
  local label="$1"
  local pidfile="$2"
  local logfile="$3"
  local pattern="$4"
  local watcher_pattern
  watcher_pattern="$(watcher_pattern_for_target "${label}")"
  shift 4
  ensure_watch_dirs
  if [ -f "${pidfile}" ]; then
    local pid
    pid="$(cat "${pidfile}")"
    if pid_running "${pid}"; then
      printf '%s already running pid=%s\n' "${label}" "${pid}"
      return
    fi
  fi
  local pid
  pid="$(find_running_pid_by_patterns "${pattern}" "${watcher_pattern}" || true)"
  if pid_running "${pid}"; then
    echo "${pid}" > "${pidfile}"
    printf '%s already running pid=%s\n' "${label}" "${pid}"
    return
  fi
  setsid "${ROOT}/benchmarks/ab_bind/run_persistent_watch.sh" "${logfile}" "$@" < /dev/null &
  local pid=$!
  echo "${pid}" > "${pidfile}"
  printf 'started %s pid=%s\n' "${label}" "${pid}"
}

stop_by_pattern() {
  local pattern="$1"
  pkill -f "${pattern}" 2>/dev/null || true
}

stop_one() {
  local label="$1"
  local pidfile="$2"
  local pattern="$3"
  local watcher_pattern
  watcher_pattern="$(watcher_pattern_for_target "${label}")"
  if [ -f "${pidfile}" ]; then
    local pid
    pid="$(cat "${pidfile}")"
    if pid_running "${pid}"; then
      kill "${pid}" 2>/dev/null || true
      printf 'stopped %s pid=%s\n' "${label}" "${pid}"
    fi
    rm -f "${pidfile}"
  fi
  stop_by_pattern "${pattern}"
  stop_by_pattern "${watcher_pattern}"
}

target_args() {
  local default_scope="$1"
  shift
  if [ "$#" -eq 0 ]; then
    case "${default_scope}" in
      start|restart)
        printf '%s\n' robust
        ;;
      stop|status)
        printf '%s\n' robust priority rescue targeted lambda sampling-qc deep ultra stale gap backlog
        ;;
    esac
    return
  fi

  local target
  for target in "$@"; do
    case "${target}" in
      robust|priority|rescue|targeted|lambda|sampling-qc|deep|ultra|stale|gap|backlog)
        printf '%s\n' "${target}"
        ;;
      all)
        printf '%s\n' robust priority rescue targeted lambda sampling-qc deep ultra stale gap backlog
        ;;
      *)
        echo "Unknown watcher target: ${target}" >&2
        echo "Valid targets: robust priority rescue targeted lambda sampling-qc deep ultra stale gap backlog all" >&2
        return 1
        ;;
    esac
  done
}

start_target() {
  local target="$1"
  case "${target}" in
    robust)
      start_one "robust" "${ROBUST_PIDFILE}" "${ROBUST_LOG}" "${ROBUST_PATTERN}" "${ROOT}/benchmarks/ab_bind/run_validation_robust_watcher.sh"
      ;;
    priority)
      start_one "priority" "${PRIORITY_PIDFILE}" "${PRIORITY_LOG}" "${PRIORITY_PATTERN}" env ABAG_ENABLE_PRIORITY_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_priority_watcher.sh"
      ;;
    rescue)
      start_one "rescue" "${RESCUE_PIDFILE}" "${RESCUE_LOG}" "${RESCUE_PATTERN}" env ABAG_ENABLE_RESCUE_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_rescue_watcher.sh"
      ;;
    targeted)
      start_one "targeted" "${TARGETED_REPEAT_SPREAD_RESCUE_PIDFILE}" "${TARGETED_REPEAT_SPREAD_RESCUE_LOG}" "${TARGETED_REPEAT_SPREAD_RESCUE_PATTERN}" env ABAG_ENABLE_TARGETED_REPEAT_SPREAD_RESCUE_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_targeted_repeat_spread_rescue_watcher.sh"
      ;;
    lambda)
      start_one "lambda" "${TARGETED_LAMBDA_RESCUE_PIDFILE}" "${TARGETED_LAMBDA_RESCUE_LOG}" "${TARGETED_LAMBDA_RESCUE_PATTERN}" env ABAG_ENABLE_TARGETED_LAMBDA_RESCUE_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_targeted_lambda_rescue_watcher.sh"
      ;;
    sampling-qc)
      start_one "sampling-qc" "${SAMPLING_QC_RESCUE_PIDFILE}" "${SAMPLING_QC_RESCUE_LOG}" "${SAMPLING_QC_RESCUE_PATTERN}" env ABAG_ENABLE_SAMPLING_QC_RESCUE_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_sampling_qc_rescue_watcher.sh"
      ;;
    deep)
      start_one "deep" "${DEEP_RESCUE_PIDFILE}" "${DEEP_RESCUE_LOG}" "${DEEP_RESCUE_PATTERN}" env ABAG_ENABLE_DEEP_RESCUE_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_deep_rescue_watcher.sh"
      ;;
    ultra)
      start_one "ultra" "${ULTRA_RESCUE_PIDFILE}" "${ULTRA_RESCUE_LOG}" "${ULTRA_RESCUE_PATTERN}" env ABAG_ENABLE_ULTRA_RESCUE_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_ultra_rescue_watcher.sh"
      ;;
    stale)
      start_one "stale" "${STALE_RECOVERY_PIDFILE}" "${STALE_RECOVERY_LOG}" "${STALE_RECOVERY_PATTERN}" env ABAG_ENABLE_STALE_RECOVERY_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_stale_recovery_watcher.sh"
      ;;
    gap)
      start_one "gap" "${GAP_RECOVERY_PIDFILE}" "${GAP_RECOVERY_LOG}" "${GAP_RECOVERY_PATTERN}" env ABAG_ENABLE_GAP_RECOVERY_WATCHER=1 "${ROOT}/benchmarks/ab_bind/run_validation_gap_recovery_watcher.sh"
      ;;
    backlog)
      start_one "backlog" "${BACKLOG_PIDFILE}" "${BACKLOG_LOG}" "${BACKLOG_PATTERN}" "${ROOT}/benchmarks/ab_bind/run_validation_priority_backlog_watcher.sh"
      ;;
  esac
}

stop_target() {
  local target="$1"
  case "${target}" in
    robust)
      stop_one "robust" "${ROBUST_PIDFILE}" "${ROBUST_PATTERN}"
      ;;
    priority)
      stop_one "priority" "${PRIORITY_PIDFILE}" "${PRIORITY_PATTERN}"
      ;;
    rescue)
      stop_one "rescue" "${RESCUE_PIDFILE}" "${RESCUE_PATTERN}"
      ;;
    targeted)
      stop_one "targeted" "${TARGETED_REPEAT_SPREAD_RESCUE_PIDFILE}" "${TARGETED_REPEAT_SPREAD_RESCUE_PATTERN}"
      ;;
    lambda)
      stop_one "lambda" "${TARGETED_LAMBDA_RESCUE_PIDFILE}" "${TARGETED_LAMBDA_RESCUE_PATTERN}"
      ;;
    sampling-qc)
      stop_one "sampling-qc" "${SAMPLING_QC_RESCUE_PIDFILE}" "${SAMPLING_QC_RESCUE_PATTERN}"
      ;;
    deep)
      stop_one "deep" "${DEEP_RESCUE_PIDFILE}" "${DEEP_RESCUE_PATTERN}"
      ;;
    ultra)
      stop_one "ultra" "${ULTRA_RESCUE_PIDFILE}" "${ULTRA_RESCUE_PATTERN}"
      ;;
    stale)
      stop_one "stale" "${STALE_RECOVERY_PIDFILE}" "${STALE_RECOVERY_PATTERN}"
      ;;
    gap)
      stop_one "gap" "${GAP_RECOVERY_PIDFILE}" "${GAP_RECOVERY_PATTERN}"
      ;;
    backlog)
      stop_one "backlog" "${BACKLOG_PIDFILE}" "${BACKLOG_PATTERN}"
      ;;
  esac
}

show_status_target() {
  local target="$1"
  case "${target}" in
    robust)
      status_line "robust" "${ROBUST_PIDFILE}" "${ROBUST_PATTERN}"
      ;;
    priority)
      status_line "priority" "${PRIORITY_PIDFILE}" "${PRIORITY_PATTERN}"
      ;;
    rescue)
      status_line "rescue" "${RESCUE_PIDFILE}" "${RESCUE_PATTERN}"
      ;;
    targeted)
      status_line "targeted" "${TARGETED_REPEAT_SPREAD_RESCUE_PIDFILE}" "${TARGETED_REPEAT_SPREAD_RESCUE_PATTERN}"
      ;;
    lambda)
      status_line "lambda" "${TARGETED_LAMBDA_RESCUE_PIDFILE}" "${TARGETED_LAMBDA_RESCUE_PATTERN}"
      ;;
    sampling-qc)
      status_line "sampling-qc" "${SAMPLING_QC_RESCUE_PIDFILE}" "${SAMPLING_QC_RESCUE_PATTERN}"
      ;;
    deep)
      status_line "deep" "${DEEP_RESCUE_PIDFILE}" "${DEEP_RESCUE_PATTERN}"
      ;;
    ultra)
      status_line "ultra" "${ULTRA_RESCUE_PIDFILE}" "${ULTRA_RESCUE_PATTERN}"
      ;;
    stale)
      status_line "stale" "${STALE_RECOVERY_PIDFILE}" "${STALE_RECOVERY_PATTERN}"
      ;;
    gap)
      status_line "gap" "${GAP_RECOVERY_PIDFILE}" "${GAP_RECOVERY_PATTERN}"
      ;;
    backlog)
      status_line "backlog" "${BACKLOG_PIDFILE}" "${BACKLOG_PATTERN}"
      ;;
  esac
}

case "${ACTION}" in
  start)
    mapfile -t TARGETS < <(target_args start "$@")
    for target in "${TARGETS[@]}"; do
      start_target "${target}"
    done
    ;;
  stop)
    mapfile -t TARGETS < <(target_args stop "$@")
    for target in "${TARGETS[@]}"; do
      stop_target "${target}"
    done
    ;;
  restart)
    mapfile -t TARGETS < <(target_args restart "$@")
    for target in "${TARGETS[@]}"; do
      stop_target "${target}"
    done
    sleep 2
    for target in "${TARGETS[@]}"; do
      start_target "${target}"
    done
    ;;
  status)
    mapfile -t TARGETS < <(target_args status "$@")
    for target in "${TARGETS[@]}"; do
      show_status_target "${target}"
    done
    ;;
  *)
    echo "Usage: $0 {start|stop|restart|status} [robust|priority|rescue|targeted|lambda|sampling-qc|deep|ultra|stale|gap|backlog|all ...]" >&2
    exit 1
    ;;
esac
