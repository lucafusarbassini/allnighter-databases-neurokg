#!/bin/bash
# =============================================================================
# Auto-backup: Continuously extract sandbox work while it runs
# =============================================================================
# Monitors the sandbox and creates backup patches every N minutes.
# Run this in a separate terminal while the sandbox is working.
#
# Usage:
#   ./auto-backup.sh [interval_minutes]
#   ./auto-backup.sh 15  # Backup every 15 minutes (default: 30)
# =============================================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INTERVAL=${1:-30}  # Default: 30 minutes
BACKUP_DIR="./backups"
mkdir -p "$BACKUP_DIR"

echo "=== Auto-Backup Started ==="
echo "Interval: ${INTERVAL} minutes"
echo "Backup dir: ${BACKUP_DIR}"
echo ""
echo "Press Ctrl+C to stop"
echo ""

while true; do
    # Check if sandbox is running
    if ! docker ps --format '{{.Names}}' | grep -q 'allnighter-sandbox'; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sandbox not running. Waiting..."
        sleep 60
        continue
    fi

    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    PATCH_FILE="${BACKUP_DIR}/auto-backup-${TIMESTAMP}.patch"
    LOG_FILE="${BACKUP_DIR}/auto-backup-${TIMESTAMP}.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating backup..."

    # Extract current diff
    docker exec allnighter-sandbox su-exec agent bash -c '
        cd /workspace
        git add -A 2>/dev/null
        git diff --cached
    ' > "$PATCH_FILE" 2>/dev/null || {
        echo "  Failed to extract diff. Skipping."
        rm -f "$PATCH_FILE"
    }

    # Extract sandbox log
    docker exec allnighter-sandbox cat /agent-logs/sandbox.log > "$LOG_FILE" 2>/dev/null || true

    # Sync data files to host (databases)
    echo "  Syncing data files..."
    mkdir -p ../template_package/data
    docker cp allnighter-sandbox:/workspace/template_package/data/. ../template_package/data/ 2>/dev/null || true

    if [ -s "$PATCH_FILE" ]; then
        LINES=$(wc -l < "$PATCH_FILE")
        echo "  Backup saved: $PATCH_FILE (${LINES} lines)"

        # Show brief summary
        docker exec allnighter-sandbox su-exec agent bash -c '
            cd /workspace
            git diff --cached --stat | tail -5
        ' 2>/dev/null || true
    else
        echo "  No changes to backup."
        rm -f "$PATCH_FILE"
    fi

    echo "  Next backup in ${INTERVAL} minutes..."
    echo ""

    # Sleep for the interval (in minutes)
    sleep $((INTERVAL * 60))
done
