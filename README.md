This project is the updated version of the FI network construction project. The original project, which is hosted at [Reactome's FI Network Project](https://github.com/reactome/fi_network_build), is based on training a Naive Bayes Classifier. This updated version is based on training random forest. There are two functions of this project:

* Generate FI files from Reactome and other manually curated pathway databases so that they can be used to generated training and test matrix files to train random forest. The project uses these files is [fi-network-ml at reactome-idg](https://github.com/reactome-idg/fi-network-ml). 


* Integrate the predicted results from the fi-network-ml project together with extracted pathway FIs to develop a FI network database and Cytoscape files so that they can be used to update ReactomeFIViz.

Note: Some of jar files need to be installed locally. Use ant scripts in the ant folder and jar files in the install_jar to do that. Please note the license statement for jacksum.jar.

The protege-related jars in the lib folder also need to be installed locally so that they can be used by maven. Run this once per machine instead of doing it by hand:

    scripts/setup_local_jars.sh

## Automating the annual build

Most of the manual procedure in `doc/ProceduresToBuildFINetwork_RF.docx` is now scripted - see the plan/scripts below rather than hand-editing `resources/configuration.prop` or running each `FINetworkBuilder` method one at a time from Eclipse:

1. `build.params` is the one place year-specific values (and your local DB password) live - it's gitignored since it contains credentials, so first run `cp build.params.example build.params` and fill it in. Edit it, then run `python3 scripts/render_config.py` to regenerate `resources/configuration.prop`, `resources/TREDHibernate.cfg.xml`, and `resources/funcIntHibernate.cfg.xml`.
2. `scripts/fetch_datasets.sh <subcommand>` (or `all`) downloads/stages the external datasets into the layout `configuration.prop` expects; `scripts/fetch_datasets.sh check-latest` reports newer upstream releases without forcing a bump. `all` only runs the datasets that actually change year to year - it skips the frozen/archival ones (`nci-pid`, `tred`, `encode`, `gene-exp-static`, `bioplex`) unless you pass `all --include-static` or run one of them by name directly.
3. `scripts/bootstrap_reactome_db.sh` creates/loads the local Reactome MySQL database and stops only at the one step that has no headless equivalent (Curator Tool's schema export).
4. `scripts/run_pipeline.sh` runs the `FINetworkBuilder` sequence headlessly in three stages, stopping between them for the two steps that need a human first (`scripts/run_pipeline.sh <stepName>` re-runs a single step directly, bypassing the staging):
   - `--stage1`: `prepareMappingFiles`, `convertPathwayDBs` - then eyeball the converted pathway projects.
   - `--stage2`: `dumpPathwayDBs`, `dumpPathwayFIs` - then run the RF training (`scripts/run_ml_local.sh`) and set `CUT_OFF_VALUE`.
   - `--stage3`: `buildFIDb`, `generateCytoscapePlugInFiles` - then runs `scripts/sanity_check.py` automatically to flag any metric that moved a lot from last year.

   All logging goes to `logs/build_<timestamp>.log`.
5. `scripts/run_ml_local.sh` wraps the local (Mac-side) RF training + precision/recall plot generation from the sibling `fi-network-ml` repo (used between stage 2 and stage 3 above).
