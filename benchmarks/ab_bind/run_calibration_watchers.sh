#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
QUICK_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_quick_plan"
RESCUE_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_calibration_rescues"
QUICK_LOG="${QUICK_ROOT}/reports/watch/watch_calibration_quick_auto.log"
RESCUE_LOG="${RESCUE_ROOT}/reports/watch/watch_calibration_rescue_auto.log"
QUICK_PIDFILE="${QUICK_ROOT}/reports/watch/watch_calibration_quick.pid"
RESCUE_PIDFILE="${RESCUE_ROOT}/reports/watch/watch_calibration_rescue.pid"

ACTION="${1:-status}"
if [ "$#" -gt 0 ]; then
  shift
fi

ensure_watch_dirs() {
  mkdir -p "$(dirname "${QUICK_LOG}")" "$(dirname "${RESCUE_LOG}")"
}

pid_running() {
  local pid="$1"
  [ -n "${pid}" ] && kill -0 "${pid}" 2>/dev/null
}

status_line() {
  local label="$1"
  local pidfile="$2"
  local pattern="$3"
  if [ -f "${pidfile}" ]; then
    local pid
    pid="$(cat "${pidfile}")"
    if pid_running "${pid}"; then
      printf '%s\tpid=%s\tstatus=running\n' "${label}" "${pid}"
      return
    fi
  fi
  local pid
  pid="$(pgrep -f "${pattern}" | head -n 1 || true)"
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
  shift 3
  ensure_watch_dirs
  if [ -f "${pidfile}" ]; then
    local pid
    pid="$(cat "${pidfile}")"
    if pid_running "${pid}"; then
      printf '%s already running pid=%s\n' "${label}" "${pid}"
      return
    fi
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
}

target_args() {
  local default_scope="$1"
  shift
  if [ "$#" -eq 0 ]; then
    case "${default_scope}" in
      start|restart)
        printf '%s\n' quick
        ;;
      stop|status)
        printf '%s\n' quick rescue
        ;;
    esac
    return
  fi

  local target
  for target in "$@"; do
    case "${target}" in
      quick|rescue)
        printf '%s\n' "${target}"
        ;;
      all)
        printf '%s\n' quick rescue
        ;;
      *)
        echo "Unknown watcher target: ${target}" >&2
        echo "Valid targets: quick rescue all" >&2
        return 1
        ;;
    esac
  done
}

start_target() {
  local target="$1"
  case "${target}" in
    quick)
      start_one \
        "calibration-quick" \
        "${QUICK_PIDFILE}" \
        "${QUICK_LOG}" \
        env ABAG_ENABLE_QUICK_WATCHER=1 \
        "${ROOT}/benchmarks/ab_bind/run_calibration_quick_watcher.sh"
      ;;
    rescue)
      start_one \
        "calibration-rescue" \
        "${RESCUE_PIDFILE}" \
        "${RESCUE_LOG}" \
        env ABAG_ENABLE_CALIBRATION_RESCUE_WATCHER=1 \
        "${ROOT}/benchmarks/ab_bind/run_calibration_rescue_watcher.sh"
      ;;
  esac
}

stop_target() {
  local target="$1"
  case "${target}" in
    quick)
      stop_one \
        "calibration-quick" \
        "${QUICK_PIDFILE}" \
        "watch_validation_priority.py --plan-root ${QUICK_ROOT} --split-name calibration"
      ;;
    rescue)
      stop_one \
        "calibration-rescue" \
        "${RESCUE_PIDFILE}" \
        "watch_validation_priority.py --plan-root ${RESCUE_ROOT} --split-name calibration"
      ;;
  esac
}

show_status_target() {
  local target="$1"
  case "${target}" in
    quick)
      status_line \
        "calibration-quick" \
        "${QUICK_PIDFILE}" \
        "watch_validation_priority.py --plan-root ${QUICK_ROOT} --split-name calibration"
      ;;
    rescue)
      status_line \
        "calibration-rescue" \
        "${RESCUE_PIDFILE}" \
        "watch_validation_priority.py --plan-root ${RESCUE_ROOT} --split-name calibration"
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
    echo "Usage: $0 {start|stop|restart|status} [quick|rescue|all ...]" >&2
    exit 1
    ;;
esac
