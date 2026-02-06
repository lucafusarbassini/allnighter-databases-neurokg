#!/bin/bash
# =============================================================================
# Watchdog: Auto-shutdown after SANDBOX_TIMEOUT_HOURS
# =============================================================================
# Provides safety net so the container doesn't run forever.
# Logs warnings at 75% and 90% of timeout before shutting down.

LOGFILE="${AGENT_LOGS}/sandbox.log"
TIMEOUT_SECONDS=$((SANDBOX_TIMEOUT_HOURS * 3600))
WARN_75=$((TIMEOUT_SECONDS * 75 / 100))
WARN_90=$((TIMEOUT_SECONDS * 90 / 100))
START_TIME=$(date +%s)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WATCHDOG] $*" | tee -a "$LOGFILE"
}

warned_75=false
warned_90=false

while true; do
    sleep 60
    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TIME))

    if [ "$ELAPSED" -ge "$WARN_75" ] && [ "$warned_75" = false ]; then
        REMAINING=$(( (TIMEOUT_SECONDS - ELAPSED) / 60 ))
        log "WARNING: 75% of timeout reached. ~${REMAINING} minutes remaining."
        warned_75=true
    fi

    if [ "$ELAPSED" -ge "$WARN_90" ] && [ "$warned_90" = false ]; then
        REMAINING=$(( (TIMEOUT_SECONDS - ELAPSED) / 60 ))
        log "WARNING: 90% of timeout reached. ~${REMAINING} minutes remaining."
        warned_90=true
    fi

    if [ "$ELAPSED" -ge "$TIMEOUT_SECONDS" ]; then
        log "TIMEOUT REACHED (${SANDBOX_TIMEOUT_HOURS}h). Shutting down sandbox."

        # Save any uncommitted work before shutdown
        if [ -d "${WORKSPACE}/.git" ]; then
            log "Saving uncommitted work..."
            cd "${WORKSPACE}" || true
            git add -A 2>/dev/null
            git commit -m "Auto-save: watchdog timeout shutdown at $(date)" 2>/dev/null || true
            log "Work saved."
        fi

        # Signal the main process to stop
        kill 1 2>/dev/null
        exit 0
    fi
done
