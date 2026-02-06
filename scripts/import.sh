#!/bin/bash
set -e  # Exit immediately if any command fails

echo "--- [IMPORT] Starting Import Script ---"

# 1. Wait for system stability
sleep 2

# 2. Import Logic
IMPORT_SCRIPT="/data/build2neo/neo4j-admin-import-call.sh"

echo "--- [IMPORT] Checking for import script at: $IMPORT_SCRIPT ---"

if [ -f "$IMPORT_SCRIPT" ]; then
    echo "--- [IMPORT] Script found. Executing... ---"
    chmod +x "$IMPORT_SCRIPT"
    "$IMPORT_SCRIPT"
else
    echo "--- [IMPORT] No import script found. Skipping import step. ---"
fi

# 3. Neo4j Lifecycle
# Your logic: Start DB, wait, Stop DB.
# This implies this container is ONLY used to initialize data, then it dies.
echo "--- [IMPORT] Starting Neo4j... ---"
neo4j start

echo "--- [IMPORT] Waiting 10 seconds for Neo4j startup/processing... ---"
sleep 10

echo "--- [IMPORT] Stopping Neo4j... ---"
neo4j stop

echo "--- [IMPORT] Script Finished. Container will now exit. ---"