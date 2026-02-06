#!/bin/bash
# =============================================================================
# Gracefully stop the sandbox
# =============================================================================
# Saves any uncommitted work before shutting down.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Stopping Sandbox ==="

# Check if running
if ! docker ps --format '{{.Names}}' | grep -q 'allnighter-sandbox'; then
    echo "Sandbox is not running."
    exit 0
fi

# Save work before stopping
echo "Saving any uncommitted work in the sandbox..."
docker exec allnighter-sandbox su-exec agent bash -c '
    cd /workspace 2>/dev/null
    git add -A 2>/dev/null
    git commit -m "Auto-save: manual shutdown at $(date)" 2>/dev/null || echo "Nothing to commit."
' || true

echo "Stopping container..."
docker compose -f docker-compose.sandbox.yml down

echo ""
echo "Sandbox stopped."
echo "To extract the work done, run: ./extract-work.sh"
echo "  (Note: you may need to restart the container first to access the volume)"
