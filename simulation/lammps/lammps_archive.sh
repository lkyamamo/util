#!/usr/bin/env bash
# =============================================================================
# lammps_archive.sh
# =============================================================================
# Copies a LAMMPS structure file (.data or .restart) to a destination
# directory and appends the lammps_describe.py output to the README.md
# "History" section in that directory.
#
# Usage:
#   lammps_archive.sh <structure_file> <lammps_input> <dest_dir> [label]
#
# Arguments:
#   structure_file  Path to the .data or .restart file to archive
#   lammps_input    Path to the LAMMPS input script (used for description)
#   dest_dir        Destination directory (created if it doesn't exist)
#   label           Optional short label/note for this entry (default: timestamp)
#
# Dependencies:
#   - lammps_describe.py  (looked up via LAMMPS_DESCRIBE env var, or alongside
#                          this script, or on $PATH)
#   - python3
#
# Examples:
#   ./lammps_archive.sh final.data in.lammps ./runs/nvt_300K "NVT equilibration at 300 K"
#   LAMMPS_DESCRIBE=/tools/lammps_describe.py ./lammps_archive.sh final.restart in.npt ./archive
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Colour helpers (suppressed if not a terminal)
# ---------------------------------------------------------------------------
if [ -t 1 ]; then
    RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
    CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; CYAN=''; BOLD=''; RESET=''
fi

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
die()     { echo -e "${RED}[ERROR]${RESET} $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    sed -n '/^# Usage/,/^# Dependencies/p' "$0" | grep '^#' | sed 's/^# \?//'
    exit 1
}

[[ $# -lt 3 ]] && usage

STRUCT_FILE="$1"
LAMMPS_INPUT="$2"
DEST_DIR="$3"
LABEL="${4:-}"

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
[[ -f "$STRUCT_FILE" ]] || die "Structure file not found: $STRUCT_FILE"
[[ -f "$LAMMPS_INPUT" ]] || die "LAMMPS input script not found: $LAMMPS_INPUT"

EXT="${STRUCT_FILE##*.}"
[[ "$EXT" == "data" || "$EXT" == "restart" ]] \
    || warn "File extension '.$EXT' is not .data or .restart — proceeding anyway."

# ---------------------------------------------------------------------------
# Locate lammps_describe.py
# ---------------------------------------------------------------------------
find_describer() {
    # 1) Explicit env var
    if [[ -n "${LAMMPS_DESCRIBE:-}" ]]; then
        [[ -f "$LAMMPS_DESCRIBE" ]] && echo "$LAMMPS_DESCRIBE" && return
        warn "LAMMPS_DESCRIBE set but not found: $LAMMPS_DESCRIBE"
    fi
    # 2) Same directory as this script
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    [[ -f "$SCRIPT_DIR/lammps_describe.py" ]] && echo "$SCRIPT_DIR/lammps_describe.py" && return
    # 3) PATH
    if command -v lammps_describe.py &>/dev/null; then
        echo "$(command -v lammps_describe.py)" && return
    fi
    echo ""
}

DESCRIBER="$(find_describer)"
if [[ -z "$DESCRIBER" ]]; then
    warn "lammps_describe.py not found — description will be skipped."
    DESCRIBER=""
fi

# ---------------------------------------------------------------------------
# Create destination directory
# ---------------------------------------------------------------------------
mkdir -p "$DEST_DIR"
success "Destination directory ready: $DEST_DIR"

# ---------------------------------------------------------------------------
# Copy structure file (never silently overwrite — version if needed)
# ---------------------------------------------------------------------------
DEST_STRUCT="$DEST_DIR/$(basename "$STRUCT_FILE")"
if [[ -e "$DEST_STRUCT" ]]; then
    STAMP="$(date +%Y%m%d_%H%M%S)"
    BASE="${STRUCT_FILE##*/}"
    NAME="${BASE%.*}"
    DEST_STRUCT="$DEST_DIR/${NAME}_${STAMP}.${EXT}"
    warn "Target exists — saving as $(basename "$DEST_STRUCT")"
fi
cp "$STRUCT_FILE" "$DEST_STRUCT"
success "Copied: $(basename "$STRUCT_FILE") → $DEST_STRUCT"

# ---------------------------------------------------------------------------
# Generate description (to a temp file)
# ---------------------------------------------------------------------------
TMP_DESC=""
if [[ -n "$DESCRIBER" ]]; then
    TMP_DESC="$(mktemp /tmp/lammps_desc_XXXXXX.txt)"
    if python3 "$DESCRIBER" "$LAMMPS_INPUT" "$TMP_DESC" 2>/dev/null; then
        success "Generated description from: $LAMMPS_INPUT"
    else
        warn "lammps_describe.py encountered an error — description may be incomplete."
    fi
fi

# ---------------------------------------------------------------------------
# Build the history entry block
# ---------------------------------------------------------------------------
TIMESTAMP="$(date '+%Y-%m-%d %H:%M:%S')"
HOSTNAME="$(hostname 2>/dev/null || echo 'unknown')"
ABS_SRC="$(cd "$(dirname "$STRUCT_FILE")" && pwd)/$(basename "$STRUCT_FILE")"
ABS_INPUT="$(cd "$(dirname "$LAMMPS_INPUT")" && pwd)/$(basename "$LAMMPS_INPUT")"
[[ -z "$LABEL" ]] && LABEL="Archived on $TIMESTAMP"

ENTRY_HEADER="### $LABEL"

build_entry() {
    echo "$ENTRY_HEADER"
    echo ""
    echo "| Field          | Value |"
    echo "|----------------|-------|"
    echo "| **Timestamp**  | \`$TIMESTAMP\` |"
    echo "| **Host**       | \`$HOSTNAME\` |"
    echo "| **Structure**  | \`$(basename "$DEST_STRUCT")\` |"
    echo "| **Source**     | \`$ABS_SRC\` |"
    echo "| **Input script** | \`$ABS_INPUT\` |"
    echo ""
    if [[ -n "$TMP_DESC" && -f "$TMP_DESC" ]]; then
        echo "<details>"
        echo "<summary>LAMMPS Input Description</summary>"
        echo ""
        echo '```'
        cat "$TMP_DESC"
        echo '```'
        echo ""
        echo "</details>"
    else
        echo "_No description available (lammps_describe.py not found)._"
    fi
    echo ""
    echo "---"
    echo ""
}

# ---------------------------------------------------------------------------
# Update README.md
# ---------------------------------------------------------------------------
README="$DEST_DIR/README.md"

# ---- Case 1: README does not exist → create it from scratch ----
if [[ ! -f "$README" ]]; then
    info "README.md not found — creating."
    {
        echo "# LAMMPS Archive"
        echo ""
        echo "This directory contains archived LAMMPS structure files and their"
        echo "associated simulation history."
        echo ""
        echo "## History"
        echo ""
        build_entry
    } > "$README"
    success "Created README.md with first history entry."

# ---- Case 2: README exists, has a ## History section ----
elif grep -q "^## History" "$README"; then
    info "Appending to existing ## History section."
    # Insert the new entry immediately after the "## History" line
    # We write to a temp file then atomically replace
    TMP_README="$(mktemp)"
    awk -v entry="$(build_entry)" '
        /^## History/ {
            print          # print the "## History" heading
            print ""       # blank line
            print entry    # new entry goes first (most-recent-first)
            next
        }
        { print }
    ' "$README" > "$TMP_README"
    mv "$TMP_README" "$README"
    success "Prepended new entry to ## History section."

# ---- Case 3: README exists, no ## History section → append the section ----
else
    info "README.md exists but has no ## History section — appending."
    {
        echo ""
        echo "## History"
        echo ""
        build_entry
    } >> "$README"
    success "Appended ## History section to existing README.md."
fi

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
[[ -n "$TMP_DESC" && -f "$TMP_DESC" ]] && rm -f "$TMP_DESC"

echo ""
echo -e "${BOLD}Done.${RESET}"
echo -e "  Structure : $DEST_STRUCT"
echo -e "  README    : $README"