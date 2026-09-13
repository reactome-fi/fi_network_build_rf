#!/usr/bin/env bash
# Wraps the local steps of the random-forest ML pipeline that live in the
# sibling reactome-idg/fi-network-ml repo: training + precision/recall/ROC plot
# generation. The remote high-memory feature-generation run (OHSU box) and the
# final CUT_OFF_VALUE pick from the plot stay manual by design - see the plan.
#
# NOTE (FI_2026): this machine's 18GB RAM is not enough to run rf_train_predict.py
# against these feature files (train ~11.7M rows / 2.6GB, prediction ~22.9M rows /
# 5.1GB) - training needs to run on a bigger machine (e.g. MacPro, per the manual).
# Use --plots-only here after copying just the (much smaller) precision_recall_*.csv
# / aggregated_test_predict_scores_*.csv / train_test_*.csv outputs back from
# wherever training actually ran.
#
# Usage:
#   scripts/run_ml_local.sh <train_feature_file> <test_feature_file> <prediction_feature_file> [down_sample] [train_test_split]
#   scripts/run_ml_local.sh --plots-only
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

if [ ! -f .build_env.sh ]; then
    echo "Missing .build_env.sh - run: python3 scripts/render_config.py" >&2
    exit 1
fi
source .build_env.sh

if [ -z "${IDG_FI_NETWORK_ML_DIR:-}" ] || [ ! -d "$IDG_FI_NETWORK_ML_DIR" ]; then
    echo "Set IDG_FI_NETWORK_ML_DIR in build.params to a local checkout of" >&2
    echo "https://github.com/reactome-idg/fi-network-ml, then re-run: python3 scripts/render_config.py" >&2
    exit 1
fi

# Prefer the conda env this repo's scripts are verified to work with
# (pandas 2.2.1 / scikit-learn 1.4.1 confirmed present) - override with
# ML_CONDA_ENV=<name> if a different one should be used.
CONDA_ENV="${ML_CONDA_ENV:-scanpy_3_10}"
PYTHON="python3"
# Check the env directory directly instead of piping `conda env list` into
# grep -q - grep exiting early on a match sends conda's process a SIGPIPE,
# which conda's own Python code doesn't handle cleanly (crashes with
# BrokenPipeError), silently defeating this whole check.
CONDA_ENV_PYTHON="$HOME/miniforge3/envs/$CONDA_ENV/bin/python"
if [ -x "$CONDA_ENV_PYTHON" ]; then
    PYTHON="$CONDA_ENV_PYTHON"
else
    echo "WARNING: conda env '$CONDA_ENV' not found at $CONDA_ENV_PYTHON, falling back to plain python3 - verify it has pandas/scikit-learn." >&2
fi

mkdir -p "$RESULT_DIR"
cd "$IDG_FI_NETWORK_ML_DIR/scripts/ml"

# F1Analyzer.analyze_auc() calls plt.show() unconditionally after saving the ROC
# png - with an interactive matplotlib backend that opens a real GUI window and
# blocks forever when run headlessly/in the background (confirmed: it hung with
# near-zero CPU use, stuck waiting on a window that can never be closed here).
# Forcing the non-interactive Agg backend makes plt.show() a no-op.
export MPLBACKEND=Agg
# plot_rf_performance() also calls plotly's fig.show() right before the
# fig.write_html() call that actually produces the output we want. BROWSER=true
# alone did not stop it from blocking (confirmed: still hung, near-zero CPU,
# after the matplotlib fix took effect) - plotly's default "browser" renderer
# doesn't go through Python's webbrowser module the way that fix assumed.
# Force the "json" renderer instead, which does no display/IO of any kind.
export PLOTLY_RENDERER=json
export BROWSER=true

if [ "${1:-}" = "--plots-only" ]; then
    echo "[run_ml_local] --plots-only: skipping training, generating plots from existing outputs in \$RESULT_DIR"
else
    TRAIN_FILE="${1:?usage: run_ml_local.sh <train_feature_file> <test_feature_file> <prediction_feature_file> [down_sample] [train_test_split]}"
    TEST_FILE="${2:?}"
    PREDICTION_FILE="${3:?}"
    DOWN_SAMPLE="${4:-none}"
    TRAIN_TEST_SPLIT="${5:-0.10}"

    # Real signature (verified against the live script, not just the manual):
    #   rf_train_predict.py <train> <test> <prediction> <working_dir> <postfix> <down_sample> <train_test_split> [no_plot]
    # All 7 positional args are required (the script itself exits if fewer are given).
    echo "[run_ml_local] Training random forest (postfix=$RF_DATE, down_sample=$DOWN_SAMPLE, train_test_split=$TRAIN_TEST_SPLIT)..."
    "$PYTHON" rf_train_predict.py "$TRAIN_FILE" "$TEST_FILE" "$PREDICTION_FILE" "$RESULT_DIR" "$RF_DATE" "$DOWN_SAMPLE" "$TRAIN_TEST_SPLIT" \
        > "$RESULT_DIR/rf_train_predict_${RF_DATE}.log" 2>&1
fi

# Main test set: F1 + ROC/AUC (needs the original test feature file for AUC).
TEST_FILE_FOR_PLOT="${2:-$RESULT_DIR/test_feature_file_${FEATURE_FILE_DATE}.csv}"
echo "[run_ml_local] Generating precision/recall + ROC plot (main test set)..."
"$PYTHON" F1Analyzer.py "$RESULT_DIR/precision_recall_${RF_DATE}.csv" "$TEST_FILE_FOR_PLOT" \
    "$RESULT_DIR/aggregated_test_predict_scores_${RF_DATE}.csv"

# 10%-holdout Reactome FI test set: F1 only, no ROC (per the manual - the feature
# file for this subset isn't separately exported).
echo "[run_ml_local] Generating F1 plot (10%-holdout Reactome FI test set)..."
"$PYTHON" F1Analyzer.py "$RESULT_DIR/train_test_precision_recall_${RF_DATE}.csv"

cat <<EOF

[run_ml_local] Plots generated under $RESULT_DIR.
Current CUT_OFF_VALUE in build.params: ${CUT_OFF_VALUE:-<blank>}
Review both plots, then set CUT_OFF_VALUE= in build.params, re-run
"python3 scripts/render_config.py", and only then run
"scripts/run_pipeline.sh --stage3".
EOF
