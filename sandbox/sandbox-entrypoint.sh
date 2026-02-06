#!/bin/bash
set -e

LOGFILE="${AGENT_LOGS}/sandbox.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

log "=== Sandbox Starting ==="
log "Timeout: ${SANDBOX_TIMEOUT_HOURS} hours"
log "Workspace: ${WORKSPACE}"

# -------------------------------------------------------
# 1. Start Docker-in-Docker daemon (runs as root via dind)
# -------------------------------------------------------
log "Starting Docker daemon (Docker-in-Docker)..."
dockerd-entrypoint.sh dockerd &>/dev/null &
DOCKERD_PID=$!

# Wait for Docker daemon to be ready
for i in $(seq 1 30); do
    if docker info &>/dev/null; then
        log "Docker daemon ready."
        break
    fi
    if [ "$i" -eq 30 ]; then
        log "WARNING: Docker daemon did not start within 30s. Pipeline commands may fail."
    fi
    sleep 1
done

# -------------------------------------------------------
# 2. Copy project files into workspace (isolation layer)
# -------------------------------------------------------
if [ -d "/project-source" ] && [ "$(ls -A /project-source 2>/dev/null)" ]; then
    log "Copying project files into workspace..."
    cp -a /project-source/. "${WORKSPACE}/"
    chown -R agent:agent "${WORKSPACE}"
    log "Project files copied."
else
    log "No project source mounted at /project-source. Workspace is empty."
fi

# -------------------------------------------------------
# 3. Copy Claude credentials with proper permissions
# -------------------------------------------------------
if [ -d "/claude-source" ] && [ "$(ls -A /claude-source 2>/dev/null)" ]; then
    log "Copying Claude credentials to agent home..."
    mkdir -p /home/agent/.claude
    cp -a /claude-source/. /home/agent/.claude/
    chown -R agent:agent /home/agent/.claude
    log "Claude credentials copied."
else
    log "WARNING: No Claude credentials mounted. Agent may not be able to authenticate."
fi

# -------------------------------------------------------
# 4. Configure git and initialize workspace repo
# -------------------------------------------------------
# Set up git identity for the agent
su-exec agent git config --global user.email "sandbox-agent@allnighter"
su-exec agent git config --global user.name "Sandbox Agent"

if [ ! -d "${WORKSPACE}/.git" ]; then
    su-exec agent git -C "${WORKSPACE}" init
    su-exec agent git -C "${WORKSPACE}" add -A
    su-exec agent git -C "${WORKSPACE}" commit -m "Sandbox baseline snapshot" 2>/dev/null || true
    log "Git repo initialized in workspace with baseline snapshot."
else
    log "Existing git repo found in workspace."
fi

# -------------------------------------------------------
# 5. Start the watchdog timer
# -------------------------------------------------------
log "Starting watchdog (auto-shutdown in ${SANDBOX_TIMEOUT_HOURS}h)..."
/usr/local/bin/watchdog.sh &
WATCHDOG_PID=$!

# -------------------------------------------------------
# 6. Drop to agent user and run Claude or keep alive
# -------------------------------------------------------
log "=== Sandbox Ready ==="
log "To attach: docker exec -it allnighter-sandbox bash"
log "To run Claude: docker exec -it allnighter-sandbox su-exec agent claude"
log ""

# If a command was passed, run it as the agent user
if [ $# -gt 0 ]; then
    log "Running command: $*"
    exec su-exec agent "$@"
else
    # Keep the container alive, waiting for interactive attachment
    log "Sandbox idle - waiting for interactive session..."
    # Trap SIGTERM for clean shutdown
    trap 'log "Received shutdown signal."; kill $DOCKERD_PID $WATCHDOG_PID 2>/dev/null; exit 0' SIGTERM SIGINT
    while true; do
        sleep 60
    done
fi
