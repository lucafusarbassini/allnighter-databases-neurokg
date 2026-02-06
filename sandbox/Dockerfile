# =============================================================================
# Claude Agent Overnight Sandbox
# =============================================================================
# Isolated environment where a Claude agent can safely operate on this
# project without risking the host system or important data.
#
# Features:
#   - Python 3.12 + Poetry + Node.js (for Claude CLI)
#   - Docker-in-Docker (agent can run the KG pipeline inside the sandbox)
#   - Non-root user with limited privileges
#   - Auto-timeout after configurable duration
#   - Activity logging
# =============================================================================

FROM docker:27-dind

# Install system dependencies
RUN apk add --no-cache \
    bash \
    curl \
    git \
    python3 \
    py3-pip \
    py3-virtualenv \
    nodejs \
    npm \
    shadow \
    su-exec \
    coreutils \
    tini

# Install Poetry
RUN pip3 install --break-system-packages poetry

# Install Claude CLI globally
RUN npm install -g @anthropic-ai/claude-code

# Create a non-root user for the agent
RUN addgroup -g 1000 agent && \
    adduser -D -u 1000 -G agent -s /bin/bash agent && \
    addgroup agent docker 2>/dev/null || true

# Create working directories
RUN mkdir -p /workspace /agent-logs && \
    chown -R agent:agent /workspace /agent-logs

# Copy entrypoint and helper scripts
COPY sandbox-entrypoint.sh /usr/local/bin/sandbox-entrypoint.sh
COPY watchdog.sh /usr/local/bin/watchdog.sh
RUN chmod +x /usr/local/bin/sandbox-entrypoint.sh /usr/local/bin/watchdog.sh

# Default environment
ENV SANDBOX_TIMEOUT_HOURS=8
ENV WORKSPACE=/workspace
ENV AGENT_LOGS=/agent-logs

# Use tini as init system for proper signal handling
ENTRYPOINT ["tini", "--", "/usr/local/bin/sandbox-entrypoint.sh"]
