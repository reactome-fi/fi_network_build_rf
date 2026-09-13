#!/usr/bin/env bash
# Thin wrapper around org.reactome.fi.PipelineRunner. Owns the three build stages
# and the two manual checkpoints between them, so the pipeline can't accidentally
# skip past a step that needs a human first:
#
#   --stage1   prepareMappingFiles, convertPathwayDBs
#              -> STOP: eyeball the converted pathway projects (see
#                 FINetworkBuilder.dumpPathwayDBs()'s javadoc)
#   --stage2   dumpPathwayDBs, dumpPathwayFIs
#              -> STOP: run scripts/run_ml_local.sh, review the precision/recall
#                 plot, set CUT_OFF_VALUE in build.params, re-run render_config.py
#   --stage3   buildFIDb, generateCytoscapePlugInFiles
#              -> runs scripts/sanity_check.py automatically at the end
#
# You can also run/re-run a single named step directly for spot fixes, e.g.:
#   scripts/run_pipeline.sh buildFIDb
#   scripts/run_pipeline.sh org.reactome.data.ReactomeDatabaseModifier.changeMyISAMToInnodb
# Bare names with no dot are prefixed with org.reactome.fi.FINetworkBuilder. -
# this bypasses staging/checkpoints, so only use it once you know the
# prerequisites for that step are already in place.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

if [ ! -f .build_env.sh ]; then
    echo "Missing .build_env.sh - run: python3 scripts/render_config.py" >&2
    exit 1
fi
source .build_env.sh

STAGE1_STEPS=(prepareMappingFiles convertPathwayDBs)
STAGE2_STEPS=(dumpPathwayDBs dumpPathwayFIs)
STAGE3_STEPS=(buildFIDb generateCytoscapePlugInFiles)

usage() {
    cat <<EOF
Usage: $0 --stage1 | --stage2 | --stage3 | <stepName> [<stepName> ...]

  --stage1   prepareMappingFiles, convertPathwayDBs
  --stage2   dumpPathwayDBs, dumpPathwayFIs
  --stage3   buildFIDb, generateCytoscapePlugInFiles

Run these in order; each of --stage1 and --stage2 stops afterward with
instructions for the manual step that has to happen before the next stage.
EOF
}

run_steps() {
    local steps=()
    for arg in "$@"; do
        if [[ "$arg" == *.*.* ]]; then
            steps+=("$arg")
        else
            steps+=("org.reactome.fi.FINetworkBuilder.$arg")
        fi
    done
    # `mvn compile exec:java` doesn't work on this machine (see scripts/compile.sh
    # for why) - compile directly and run with java instead.
    ./scripts/compile.sh
    local cp
    cp=$(cat target/runtime_classpath.txt)
    # Several steps (convertPathwayDBs, copyHumanReferenceGeneProducts, buildFIDb)
    # are documented as needing 10-12G heap for realistic data volumes.
    java -Xmx"${JAVA_XMX:-10G}" -cp "$cp" org.reactome.fi.PipelineRunner "${steps[@]}"
}

if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

case "$1" in
    --stage1)
        run_steps "${STAGE1_STEPS[@]}"
        cat <<EOF

===== STAGE 1 DONE - MANUAL STEP BEFORE STAGE 2 =====
Eyeball the converted pathway projects under \$DATA_SET_DIR (NCI-PID, Panther,
TRED, ENCODE .rtpj files) to make sure nothing looks weird - see
FINetworkBuilder.dumpPathwayDBs()'s javadoc.

Then run: scripts/run_pipeline.sh --stage2
======================================================
EOF
        ;;
    --stage2)
        run_steps "${STAGE2_STEPS[@]}"
        cat <<EOF

===== STAGE 2 DONE - MANUAL STEP BEFORE STAGE 3 =====
buildFIDb needs RF_PREDICTION_FILE/RF_FEATURE_FILE, which don't exist yet:
  1. Run the feature-generation jar on the remote high-memory box (out of
     scope here - see the manual), copy the resulting feature files back.
  2. Run scripts/run_ml_local.sh to train locally and generate the
     precision/recall plot.
  3. Review the plot, set CUT_OFF_VALUE in build.params, then re-run:
       python3 scripts/render_config.py

Then run: scripts/run_pipeline.sh --stage3
======================================================
EOF
        ;;
    --stage3)
        run_steps "${STAGE3_STEPS[@]}"
        python3 scripts/sanity_check.py || true
        ;;
    --help|-h)
        usage
        ;;
    *)
        run_steps "$@"
        ;;
esac
