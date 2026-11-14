#!/bin/bash
# One-time setup for a checkout: point git at the tracked hooks in .githooks/.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
chmod +x .githooks/* tools/*.sh tools/*.py
git config core.hooksPath .githooks
echo "hooks enabled: $(git config core.hooksPath)"
