# Sourced by the hooks: sets PY to the first python3 that can import sqlite3
# (the ledger tests need it) and REPO_CHECK to the gate. A pyenv build linked
# against a removed MacPorts libsqlite fails the import and is skipped.
# Set UTIL_PYTHON to choose the interpreter.

REPO_CHECK="$(git rev-parse --show-toplevel)/tools/repo_check.py"
PY=""
for candidate in "${UTIL_PYTHON:-}" python3.14 python3.13 python3.12 python3.11 python3 /usr/bin/python3; do
  [ -n "$candidate" ] || continue
  if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sqlite3' 2>/dev/null; then
    PY="$candidate"
    break
  fi
done
if [ -z "$PY" ]; then
  echo "git hook: no python3 that can import sqlite3; set UTIL_PYTHON" >&2
  exit 1
fi
