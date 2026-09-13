#!/usr/bin/env bash
# Bootstraps the local Reactome MySQL database used as the FI network build's data
# source: creates the DBs, loads the dump, converts MyISAM -> InnoDB, copies missing
# human ReferenceGeneProducts, and applies the schema modification SQL. Stops right
# before the one step that has no headless equivalent (Curator Tool's schema export).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

if [ ! -f .build_env.sh ]; then
    echo "Missing .build_env.sh - run: python3 scripts/render_config.py" >&2
    exit 1
fi
source .build_env.sh

log() { echo "[bootstrap_reactome_db] $*"; }
mysql_cmd() { mysql -u"$DB_USER" -p"$DB_PWD" "$@"; }

START_TS=$(date +%s)

log "Creating databases (if not present): $REACTOME_SOURCE_DB_NAME, $REACTOME_GK_CENTRAL_DB_NAME"
mysql_cmd -e "CREATE DATABASE IF NOT EXISTS $REACTOME_SOURCE_DB_NAME"
mysql_cmd -e "CREATE DATABASE IF NOT EXISTS $REACTOME_GK_CENTRAL_DB_NAME"

DUMP="$DATASET_ROOT/reactome_dump/gk_current.sql.gz"
if [ ! -f "$DUMP" ]; then
    echo "Missing $DUMP - run: scripts/fetch_datasets.sh reactome-dump" >&2
    exit 1
fi

TABLE_COUNT=$(mysql_cmd -N -e "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema='$REACTOME_SOURCE_DB_NAME'")
if [ "$TABLE_COUNT" -eq 0 ]; then
    log "Loading gk_current.sql.gz into $REACTOME_SOURCE_DB_NAME (this can take a while)..."
    # gunzip -c (not zcat - macOS's zcat expects .Z compress-format, not .gz, and
    # fails silently with an empty pipe when given a .gz file) works portably.
    gunzip -c "$DUMP" | mysql_cmd "$REACTOME_SOURCE_DB_NAME"
else
    log "$REACTOME_SOURCE_DB_NAME already has $TABLE_COUNT tables, skipping dump load."
fi

GK_CENTRAL_DUMP="$DATASET_ROOT/reactome_dump/gk_central.sql.gz"
if [ -f "$GK_CENTRAL_DUMP" ]; then
    GK_TABLE_COUNT=$(mysql_cmd -N -e "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema='$REACTOME_GK_CENTRAL_DB_NAME'")
    if [ "$GK_TABLE_COUNT" -eq 0 ]; then
        log "Loading gk_central snapshot into $REACTOME_GK_CENTRAL_DB_NAME..."
        gunzip -c "$GK_CENTRAL_DUMP" | mysql_cmd "$REACTOME_GK_CENTRAL_DB_NAME"
    else
        log "$REACTOME_GK_CENTRAL_DB_NAME already has $GK_TABLE_COUNT tables, skipping."
    fi
else
    log "NOTE: $GK_CENTRAL_DUMP not found. The gk_central snapshot (curator.reactome.org) requires"
    log "      curator credentials and can't be scripted here - get a dump from a maintainer, place"
    log "      it at that path, and re-run this script to load it before continuing."
fi

log "Converting MyISAM to InnoDB and dropping full-text indexes..."
./scripts/run_pipeline.sh org.reactome.data.ReactomeDatabaseModifier.changeMyISAMToInnodb

log "Copying missing human ReferenceGeneProducts from the gk_central snapshot..."
./scripts/run_pipeline.sh org.reactome.data.ReactomeDatabaseModifier.copyHumanReferenceGeneProducts

log "Applying schema modification (Interaction/TargetedInteraction classes)..."
mysql_cmd "$REACTOME_SOURCE_DB_NAME" < resources/SchemaModification.sql

cat <<EOF

=====================================================================
MANUAL STEP REQUIRED (Curator Tool has no headless schema export):
  1. Open Curator Tool
  2. Connect to: localhost:3306 / $REACTOME_SOURCE_DB_NAME / $DB_USER
  3. File > Export Schema
  4. Save as exactly: $(pwd)/resources/schema   (no file extension)
=====================================================================

EOF

while true; do
    read -r -p "Press Enter once resources/schema has been saved..." _
    if [ -f resources/schema ]; then
        MTIME=$(stat -f %m resources/schema 2>/dev/null || stat -c %Y resources/schema)
        if [ "$MTIME" -ge "$START_TS" ]; then
            break
        fi
    fi
    echo "resources/schema not found or not updated since this script started - please export it now."
done

log "Reactome DB bootstrap complete."
