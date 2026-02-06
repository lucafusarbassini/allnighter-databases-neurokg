#!/bin/bash
# =============================================================================
# Extract work from the sandbox back to the host
# =============================================================================
# Copies the agent's workspace changes as a git patch that you can review
# and selectively apply to your actual project.
#
# Usage:
#   ./extract-work.sh                  # Generate patch file
#   ./extract-work.sh --apply          # Generate and apply patch
# =============================================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
PATCH_DIR="${SCRIPT_DIR}/patches"
PATCH_FILE="${PATCH_DIR}/sandbox_work_${TIMESTAMP}.patch"
APPLY=false

if [ "$1" = "--apply" ]; then
    APPLY=true
fi

# Check if container is running
if ! docker ps --format '{{.Names}}' | grep -q 'allnighter-sandbox'; then
    echo "ERROR: Sandbox container is not running."
    echo "If it was stopped, you can still extract from the volume:"
    echo "  docker compose -f docker-compose.sandbox.yml start sandbox"
    exit 1
fi

mkdir -p "$PATCH_DIR"

echo "=== Extracting work from sandbox ==="

# Generate a diff of all changes since the baseline snapshot
echo "Generating patch..."
docker exec allnighter-sandbox su-exec agent bash -c '
    cd /workspace
    git add -A 2>/dev/null
    git diff --cached
' > "$PATCH_FILE"

if [ ! -s "$PATCH_FILE" ]; then
    echo "No changes detected in the sandbox workspace."
    rm -f "$PATCH_FILE"
    exit 0
fi

LINES=$(wc -l < "$PATCH_FILE")
echo "Patch generated: ${PATCH_FILE} (${LINES} lines)"
echo ""

# Show summary of changed files
echo "Files changed:"
docker exec allnighter-sandbox su-exec agent bash -c '
    cd /workspace
    git diff --cached --stat
'
echo ""

# Also extract the sandbox log
LOG_FILE="${PATCH_DIR}/sandbox_log_${TIMESTAMP}.log"
docker exec allnighter-sandbox cat /agent-logs/sandbox.log > "$LOG_FILE" 2>/dev/null || true
echo "Sandbox log: ${LOG_FILE}"

if [ "$APPLY" = true ]; then
    echo ""
    echo "Applying patch to project..."
    cd "$PROJECT_DIR"
    git apply --stat "$PATCH_FILE"
    echo ""
    read -p "Apply these changes? [y/N] " confirm
    if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
        git apply "$PATCH_FILE"
        echo "Changes applied successfully."
    else
        echo "Cancelled. Patch saved at: ${PATCH_FILE}"
    fi
else
    echo ""
    echo "To review the patch:"
    echo "  less ${PATCH_FILE}"
    echo ""
    echo "To apply the patch to your project:"
    echo "  cd ${PROJECT_DIR}"
    echo "  git apply ${PATCH_FILE}"
    echo ""
    echo "Or run: ./extract-work.sh --apply"
fi
