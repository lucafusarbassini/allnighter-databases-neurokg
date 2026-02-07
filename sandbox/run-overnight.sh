#!/bin/bash
# =============================================================================
# Run Claude overnight loop inside the sandbox
# =============================================================================
# Each iteration has a hard timeout (default 90 min) to prevent stalls.
# If Claude gets stuck (context window, API errors, etc.), the timeout
# kills it and the loop advances to the next iteration.
#
# Usage (from host):
#   docker exec -it allnighter-sandbox su-exec agent /workspace/sandbox/run-overnight.sh
#   docker exec -it allnighter-sandbox su-exec agent /workspace/sandbox/run-overnight.sh 20 90
#
# Args:
#   $1 = number of iterations (default: 20)
#   $2 = timeout per iteration in minutes (default: 90)
# =============================================================================

set -e
cd /workspace

ITERATIONS=${1:-20}
TIMEOUT_MIN=${2:-90}
PROMPT_FILE="/workspace/martinprompt.md"

if [ ! -f "$PROMPT_FILE" ]; then
    echo "ERROR: $PROMPT_FILE not found!"
    exit 1
fi

PROMPT=$(cat "$PROMPT_FILE")

echo "=== Overnight Loop Starting ==="
echo "Iterations: ${ITERATIONS}"
echo "Timeout per iteration: ${TIMEOUT_MIN} minutes"
echo ""

for i in $(seq 1 "$ITERATIONS"); do
    echo "=== Iteration ${i}/${ITERATIONS} ==="
    START=$(date +%s)

    # Run Claude with a hard timeout
    timeout "${TIMEOUT_MIN}m" claude --dangerously-skip-permissions --model opus -p "$PROMPT" || {
        EXIT_CODE=$?
        if [ "$EXIT_CODE" -eq 124 ]; then
            echo ""
            echo "[LOOP] Iteration ${i} TIMED OUT after ${TIMEOUT_MIN} minutes. Moving to next."
        else
            echo ""
            echo "[LOOP] Iteration ${i} exited with code ${EXIT_CODE}."
        fi
    }

    ELAPSED=$(( ($(date +%s) - START) / 60 ))
    echo "[LOOP] Iteration ${i} ran for ${ELAPSED} minutes."

    # Auto-commit any uncommitted work after each iteration
    git add -A 2>/dev/null
    git commit -m "Auto-commit after iteration ${i} (${ELAPSED}min)" 2>/dev/null || true

    # Brief pause between iterations
    sleep 5
done

echo "=== Overnight Loop Complete: ${ITERATIONS} iterations ==="
