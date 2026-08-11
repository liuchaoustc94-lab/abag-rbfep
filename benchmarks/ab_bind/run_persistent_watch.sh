#!/usr/bin/env bash
set -uo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <logfile> <command> [args...]" >&2
  exit 1
fi

LOGFILE="$1"
shift

RESTART_DELAY_SECONDS="${WATCH_RESTART_DELAY_SECONDS:-5}"
NO_RESTART_CODE="${WATCH_SUPERVISOR_NO_RESTART_CODE:-75}"
mkdir -p "$(dirname "${LOGFILE}")"

child_pid=""

cleanup() {
  if [ -n "${child_pid}" ] && kill -0 "${child_pid}" 2>/dev/null; then
    kill -- "-${child_pid}" 2>/dev/null || kill "${child_pid}" 2>/dev/null || true
    wait "${child_pid}" 2>/dev/null || true
  fi
  exit 0
}

trap cleanup TERM INT HUP

while true; do
  printf '[watch-supervisor] starting at %s: %s\n' "$(date -Iseconds)" "$*" >>"${LOGFILE}"
  setsid "$@" >>"${LOGFILE}" 2>&1 < /dev/null &
  child_pid=$!
  wait "${child_pid}"
  rc=$?
  printf '[watch-supervisor] child exited code=%s at %s\n' "${rc}" "$(date -Iseconds)" >>"${LOGFILE}"
  child_pid=""
  if [ "${rc}" -eq "${NO_RESTART_CODE}" ]; then
    printf '[watch-supervisor] not restarting after disable code=%s at %s\n' "${rc}" "$(date -Iseconds)" >>"${LOGFILE}"
    exit 0
  fi
  sleep "${RESTART_DELAY_SECONDS}"
done
