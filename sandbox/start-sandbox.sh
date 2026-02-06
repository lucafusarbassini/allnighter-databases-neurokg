#!/bin/bash
# =============================================================================
# Start the Claude Agent Sandbox
# =============================================================================
# Usage:
#   ./start-sandbox.sh              # Start with defaults (8h timeout)
#   ./start-sandbox.sh --hours 12   # Custom timeout
#   ./start-sandbox.sh --attach     # Start and immediately attach
# =============================================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Defaults
TIMEOUT_HOURS=8
ATTACH=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --hours)
            TIMEOUT_HOURS="$2"
            shift 2
            ;;
        --attach)
            ATTACH=true
            shift
            ;;
        --help)
            echo "Usage: $0 [--hours N] [--attach]"
            echo ""
            echo "Options:"
            echo "  --hours N    Set sandbox timeout in hours (default: 8)"
            echo "  --attach     Attach to Claude session immediately after start"
            echo ""
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check for .env file
if [ ! -f ".env" ]; then
    echo "ERROR: .env file not found."
    echo "Copy the example and fill in your API key:"
    echo "  cp .env.example .env"
    echo "  # Then edit .env with your ANTHROPIC_API_KEY"
    exit 1
fi

# Source .env for validation
source .env
if [ -z "$ANTHROPIC_API_KEY" ] || [ "$ANTHROPIC_API_KEY" = "sk-ant-xxxxxxxxxxxxxxxx" ]; then
    echo "ERROR: ANTHROPIC_API_KEY not set in .env"
    echo "Edit .env and add your actual API key."
    exit 1
fi

# Export timeout
export SANDBOX_TIMEOUT_HOURS="$TIMEOUT_HOURS"

echo "=== Allnighter Sandbox ==="
echo "Timeout:   ${TIMEOUT_HOURS} hours"
echo "Workspace: isolated copy of project"
echo ""

# Build and start
echo "Building sandbox image..."
docker compose -f docker-compose.sandbox.yml build

echo "Starting sandbox container..."
docker compose -f docker-compose.sandbox.yml up -d

echo ""
echo "Sandbox is running!"
echo ""
echo "Commands:"
echo "  # Attach to sandbox and run Claude interactively:"
echo "  docker exec -it allnighter-sandbox su-exec agent claude"
echo ""
echo "  # Or run Claude with a specific prompt:"
echo "  docker exec -it allnighter-sandbox su-exec agent claude -p 'your task here'"
echo ""
echo "  # Open a shell in the sandbox:"
echo "  docker exec -it allnighter-sandbox su-exec agent bash"
echo ""
echo "  # View sandbox logs:"
echo "  docker exec allnighter-sandbox cat /agent-logs/sandbox.log"
echo ""
echo "  # Extract work when done:"
echo "  ./extract-work.sh"
echo ""
echo "  # Stop the sandbox:"
echo "  docker compose -f docker-compose.sandbox.yml down"
echo ""

if [ "$ATTACH" = true ]; then
    echo "Attaching to sandbox..."
    sleep 3  # Give container time to initialize
    docker exec -it allnighter-sandbox su-exec agent claude
fi
