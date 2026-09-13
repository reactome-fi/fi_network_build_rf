#!/usr/bin/env bash
# Downloads/stages the external datasets needed by the FI network build, into the
# exact directory layout resources/configuration.prop expects (via .build_env.sh).
#
# Usage: scripts/fetch_datasets.sh <subcommand>
#   Real downloads:      uniprot pir reactome-dump pfam panther string biogrid bioplex go ensembl
#   Staged from archive:  nci-pid tred encode gene-exp-static
#   Helper:               check-latest   all
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

if [ ! -f .build_env.sh ]; then
    echo "Missing .build_env.sh - run: python3 scripts/render_config.py" >&2
    exit 1
fi
source .build_env.sh

log() { echo "[fetch_datasets] $*"; }
fail() { echo "[fetch_datasets] ERROR: $*" >&2; exit 1; }

# Downloads $1 into $2 unless $2 already exists. Aborts loudly on a 404/dead link
# instead of silently saving an HTML error page as data. Also catches the case
# curl's --fail can't: a redirect that lands on a real HTML page with a 200 status
# (e.g. a download portal's homepage) instead of the actual file.
fetch() {
    local url="$1" dest="$2"
    if [ -f "$dest" ]; then
        log "already have $dest, skipping"
        return 0
    fi
    mkdir -p "$(dirname "$dest")"
    log "downloading $url -> $dest"
    curl --fail --location --silent --show-error -o "$dest" "$url" || fail "download failed: $url"
    if head -c 256 "$dest" 2>/dev/null | grep -qi "<!doctype html\|<html"; then
        rm -f "$dest"
        fail "download landed on an HTML page instead of the expected file: $url"
    fi
}

gunzip_if_present() {
    local gz="$1"
    if [ -f "$gz" ]; then
        gunzip -kf "$gz"
    fi
}

# UniProt's FTP only ever serves "current_release" - there is no URL for a specific
# past release. So after a fresh download, record + verify what we actually got
# against the UNIPROT_RELEASE label in build.params instead of trusting the label
# blindly (the label can't be used to *request* a version, only to record one).
check_uniprot_release() {
    local dir="$1"
    local notes actual actual_label
    notes=$(curl -fsSL https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt 2>/dev/null || true)
    actual=$(echo "$notes" | grep -oE 'Release [0-9]{4}_[0-9]{2}' | head -1 | awk '{print $2}')
    if [ -z "$actual" ]; then
        log "WARNING: could not determine the actual UniProt release from relnotes.txt - verify manually."
        return 0
    fi
    actual_label="release_${actual}"
    {
        echo "actual_release=$actual_label"
        echo "downloaded_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "$dir/RELEASE_INFO.txt"
    log "UniProt current_release is actually $actual_label (recorded in $dir/RELEASE_INFO.txt)"
    if [ "$actual_label" != "$UNIPROT_RELEASE" ]; then
        log "WARNING: build.params has UNIPROT_RELEASE=$UNIPROT_RELEASE but the files just"
        log "         downloaded into $dir are actually from $actual_label. Update"
        log "         UNIPROT_RELEASE in build.params to $actual_label and re-run"
        log "         'python3 scripts/render_config.py' so the label matches what's on disk."
    fi
}

cmd_uniprot() {
    local dir="$UNIPROT_DIR"
    local base="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions"
    local sprot_gz="$dir/uniprot_sprot_human.dat.gz"
    local sprot_existed=0
    [ -f "$sprot_gz" ] && sprot_existed=1

    fetch "$base/uniprot_sprot_human.dat.gz" "$sprot_gz"
    fetch "$base/uniprot_trembl_human.dat.gz" "$dir/uniprot_trembl_human.dat.gz"
    fetch "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz" \
        "$dir/uniprot_sprot_varsplic.fasta.gz"
    gunzip_if_present "$sprot_gz"
    gunzip_if_present "$dir/uniprot_trembl_human.dat.gz"
    gunzip_if_present "$dir/uniprot_sprot_varsplic.fasta.gz"

    if [ "$sprot_existed" -eq 0 ]; then
        check_uniprot_release "$dir"
    elif [ ! -f "$dir/RELEASE_INFO.txt" ]; then
        log "NOTE: $sprot_gz already present but no $dir/RELEASE_INFO.txt exists - remove"
        log "      $sprot_gz and re-run this command to record/verify its actual release."
    fi
    log "UniProt files ready in $dir"
}

cmd_pir() {
    local dir="$DATA_SET_DIR/iproclass/$DATE"
    # PIR's FTP path is known to move over time (noted in the manual) - if this
    # 404s, check https://proteininformationresource.org/pirwww/download/ for the
    # current address and update the URL below.
    fetch "ftp://ftp.proteininformationresource.org/databases/idmapping/mapping_by_sp/h_sapiens.tb" \
        "$dir/h_sapiens.tb"
    log "PIR iproclass mapping ready in $dir"
}

cmd_reactome_dump() {
    local dir="$DATASET_ROOT/reactome_dump"
    fetch "https://reactome.org/download/current/databases/gk_current.sql.gz" "$dir/gk_current.sql.gz"
    log "Reactome dump ready at $dir/gk_current.sql.gz - load it with scripts/bootstrap_reactome_db.sh"
}

cmd_pfam() {
    local dir="$PFAM_DIR_NAME"
    # Folder naming is "Pfam<release>", not just "<release>".
    local base="https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam${PFAM_RELEASE}/database_files"
    mkdir -p "$dir"
    if [ -f "$dir/pfamA_interactions.txt" ]; then
        log "Pfam files already present in $dir, skipping."
    elif curl -fsIL "$base/pfamA_interactions.txt.gz" >/dev/null 2>&1; then
        fetch "$base/pfamA.txt.gz" "$dir/pfamA.txt.gz"
        fetch "$base/pfamA_interactions.txt.gz" "$dir/pfamA_interactions.txt.gz"
        gunzip_if_present "$dir/pfamA.txt.gz"
        gunzip_if_present "$dir/pfamA_interactions.txt.gz"
    else
        # As of this writing EBI no longer hosts the database_files dump for older
        # Pfam releases at all (confirmed: 404 for release 33.1, matching the
        # manual's own suspicion that Pfam stopped maintaining domain-domain
        # interaction data). Fall back to a previous year's already-downloaded copy.
        local prev_copy
        prev_copy=$(find "$DATASET_ROOT" -maxdepth 4 -path "*/Pfam/${PFAM_RELEASE}/pfamA_interactions.txt" 2>/dev/null | head -1)
        if [ -n "$prev_copy" ]; then
            local prev_dir
            prev_dir=$(dirname "$prev_copy")
            cp "$prev_dir/pfamA.txt" "$prev_dir/pfamA_interactions.txt" "$dir/"
            [ -f "$prev_dir/version.txt" ] && cp "$prev_dir/version.txt" "$dir/"
            log "EBI no longer hosts Pfam ${PFAM_RELEASE}'s database_files - reused the copy already at $prev_dir"
        else
            fail "Pfam ${PFAM_RELEASE}'s database_files are no longer hosted by EBI, and no" \
                 "previous local copy was found under $DATASET_ROOT/*/Pfam/${PFAM_RELEASE}/." \
                 "Obtain pfamA.txt and pfamA_interactions.txt from an archived copy and place" \
                 "them in $dir."
        fi
    fi
    # PfamAnalyzer.convertIntToPfamIDs() expects this symlink name (manual step it replaces).
    (cd "$dir" && ln -sf pfamA_interactions.txt IntPFamIDs.txt)
    log "Pfam release $PFAM_RELEASE ready in $dir"
}

cmd_panther() {
    local dir="$PANTHER_DIR"
    mkdir -p "$dir/SBML"

    # Ortholog mapping (cross-species PPI -> human) lives under its own release
    # subfolder, per-year under DATA_SET_DIR - matches the real production
    # application.properties convention (panther.orthologous.map) and last year's
    # actual local layout (FI_2025/Panther/ortholog_19.0/).
    mkdir -p "$PANTHER_ORTHOLOG_DIR"
    if [ -f "$PANTHER_ORTHOLOG_DIR/HUMAN_RefGenomeOrthologs" ]; then
        log "Panther ortholog mapping already present in $PANTHER_ORTHOLOG_DIR, skipping."
    else
        # NOTE: like UniProt, this URL always serves whatever PantherDB currently
        # calls "current_release" - there is no way to fetch a specific past
        # ortholog release by URL. If PANTHER_ORTHOLOG_VERSION in build.params
        # doesn't match what's actually live, the label will be wrong; verify
        # manually against https://data.pantherdb.org/ftp/ortholog/ if unsure.
        fetch "https://data.pantherdb.org/ftp/ortholog/current_release/RefGenomeOrthologs.tar.gz" \
            "$PANTHER_ORTHOLOG_DIR/RefGenomeOrthologs.tar.gz"
        (cd "$PANTHER_ORTHOLOG_DIR" && tar -xzf RefGenomeOrthologs.tar.gz)
        grep '^HUMAN' "$PANTHER_ORTHOLOG_DIR/RefGenomeOrthologs" > "$PANTHER_ORTHOLOG_DIR/HUMAN_RefGenomeOrthologs" \
            2>/dev/null || log "NOTE: could not find RefGenomeOrthologs after extraction - check the archive layout and adjust this command"
        log "Panther ortholog mapping (release $PANTHER_ORTHOLOG_VERSION) ready in $PANTHER_ORTHOLOG_DIR"
    fi

    # Panther pathway files (SBML/BioPAX) at ftp.pantherdb.org/pathway/current_release/ have
    # changed format across releases (BioPAX-only as of recent releases), and pathways
    # haven't changed much version to version (see build.params comment) - so if this
    # PANTHER_VERSION's SBML files were already downloaded for a previous year, reuse
    # them instead of needing a fresh manual download every time.
    if [ -n "$(ls -A "$dir/SBML" 2>/dev/null)" ]; then
        log "Panther SBML files already present in $dir/SBML, skipping."
    else
        local prev_sbml
        prev_sbml=$(find "$DATASET_ROOT" -maxdepth 4 -type d -path "*/Panther/Version${PANTHER_VERSION}/SBML" ! -path "$dir/SBML" 2>/dev/null | head -1)
        if [ -n "$prev_sbml" ] && [ -n "$(ls -A "$prev_sbml" 2>/dev/null)" ]; then
            cp "$prev_sbml"/*.xml "$dir/SBML/"
            local prev_mapping
            prev_mapping=$(dirname "$prev_sbml")/SequenceAssociationPathway${PANTHER_VERSION}.txt
            [ -f "$prev_mapping" ] && cp "$prev_mapping" "$dir/"
            log "Reused Panther $PANTHER_VERSION SBML files (and mapping file, if found) from $prev_sbml"
        else
            log "NOTE: no SBML files for Panther $PANTHER_VERSION found locally or live. Download the"
            log "      pathway archive for version $PANTHER_VERSION manually into $dir/SBML and verify"
            log "      PantherToReactomeConverterTest still parses it (see build.params comment)."
        fi
    fi
}

cmd_string() {
    local base_dir="$DATA_SET_DIR/StringDB"
    local version="$STRING_VERSION"
    local base="https://stringdb-downloads.org/download"
    # species:taxon_id:SubfolderName. Not a bash associative array on purpose -
    # the default /bin/bash on macOS is 3.2, which predates `declare -A`.
    # Subfolder-per-species with uncompressed .txt matches the layout the
    # idg-fi-network-ml project's FeatureFileGenerator actually expects (verified
    # against last year's real FI_2024/StringDB/ layout, not just the manual).
    local species_taxa="human:9606:Human yeast:4932:Yeast fly:7227:Fly worm:6239:Worm mouse:10090:Mouse"
    local pair species taxon subdir dir gz
    for pair in $species_taxa; do
        species="${pair%%:*}"
        local rest="${pair#*:}"
        taxon="${rest%%:*}"
        subdir="${rest#*:}"
        dir="$base_dir/$subdir"
        gz="$dir/${taxon}.protein.links.full.${version}.txt.gz"
        fetch "$base/protein.links.full.${version}/${taxon}.protein.links.full.${version}.txt.gz" "$gz"
        gunzip_if_present "$gz"
        if [ "$species" = "human" ]; then
            gz="$dir/${taxon}.protein.info.${version}.txt.gz"
            fetch "$base/protein.info.${version}/${taxon}.protein.info.${version}.txt.gz" "$gz"
            gunzip_if_present "$gz"
        else
            gz="$dir/${taxon}.protein.aliases.${version}.txt.gz"
            fetch "$base/protein.aliases.${version}/${taxon}.protein.aliases.${version}.txt.gz" "$gz"
            gunzip_if_present "$gz"
        fi
    done
    log "StringDB PPI files ready in $base_dir (Human/Mouse/Fly/Worm/Yeast subfolders, uncompressed)"
}

cmd_biogrid() {
    local dir="$DATA_SET_DIR/BioGrid"
    # The plain Release-Archive path now 302s to the HTML download portal instead of
    # serving the file - the real static-file path has an extra /Download/ segment
    # (confirmed by inspecting the portal page's actual download link).
    local base="https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-${BIOGRID_VERSION}"
    # Per-species org files go under BIOGRID-ORGANISM-<major>/ (e.g. BIOGRID-ORGANISM-4
    # for 4.4.232); the identifiers file stays flat. Matches the real production
    # application.properties convention (biogrid.dir/biogrid.id.file.selected).
    local major="${BIOGRID_VERSION%%.*}"
    local org_dir="$dir/BIOGRID-ORGANISM-${major}"
    mkdir -p "$org_dir"
    fetch "$base/BIOGRID-ORGANISM-${BIOGRID_VERSION}.tab2.zip" "$dir/BIOGRID-ORGANISM-${BIOGRID_VERSION}.tab2.zip"
    fetch "$base/BIOGRID-IDENTIFIERS-${BIOGRID_VERSION}.tab.zip" "$dir/BIOGRID-IDENTIFIERS-${BIOGRID_VERSION}.tab.zip"
    (cd "$org_dir" && unzip -o -q "$dir/BIOGRID-ORGANISM-${BIOGRID_VERSION}.tab2.zip")
    (cd "$dir" && unzip -o -q "BIOGRID-IDENTIFIERS-${BIOGRID_VERSION}.tab.zip")
    log "BioGRID $BIOGRID_VERSION ready in $dir (organism files under BIOGRID-ORGANISM-${major}/)"
}

cmd_bioplex() {
    local dir="$DATA_SET_DIR/BioPlex"
    fetch "https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv" \
        "$dir/BioPlex_293T_Network_10K_Dec_2019.tsv"
    fetch "https://bioplex.hms.harvard.edu/data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv" \
        "$dir/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"
    log "BioPlex files ready in $dir (these rarely change per the manual - verify before re-running)"
}

cmd_go() {
    local dir="$DATA_SET_DIR/GO/$DATE"
    fetch "https://purl.obolibrary.org/obo/go.obo" "$dir/go.obo"
    fetch "https://current.geneontology.org/annotations/goa_human.gaf.gz" "$dir/goa_human.gaf.gz"
    gunzip_if_present "$dir/goa_human.gaf.gz"
    log "GO files ready in $dir"
}

cmd_ensembl() {
    local dir="$ENSEMBL_DIR"
    mkdir -p "$dir"
    log "Ensembl compara database (ensembl_compara_${ENSEMBL_RELEASE}) is a MySQL database, not a flat"
    log "file - import it from Ensembl's public MySQL server (see ensembl.org/info/data/mysql.html) for"
    log "release $ENSEMBL_RELEASE. Place ProteinFamilies.txt at $dir/ProteinFamilies.txt."
}

cmd_nci_pid() {
    local archive="data_archive/nci_pid"
    local dir="$DATA_SET_DIR/NCI-Pathways/01162012"
    if [ -f "$archive/NCI-Nature_Curated.bp2.owl" ] && [ -f "$archive/BioCarta.bp2.owl" ]; then
        mkdir -p "$dir"
        cp -n "$archive"/*.owl "$dir/"
        log "Staged NCI-PID BioPAX files into $dir"
    else
        fail "NCI-PID BioPAX files not found in $archive/. This source has been defunct since ~2014" \
             "and cannot be re-downloaded - obtain NCI-Nature_Curated.bp2.owl and BioCarta.bp2.owl" \
             "(ask a project maintainer for an existing copy) and place them in $archive/, then re-run."
    fi
}

cmd_tred() {
    # TRED_DIR is where TREDToReactomeConverter writes its output .rtpj project file -
    # needed regardless of whether the MySQL database itself still needs loading.
    mkdir -p "$DATA_SET_DIR/TRED"
    if mysql -u"$DB_USER" -p"$DB_PWD" -e 'USE TRED' 2>/dev/null; then
        log "TRED database already loaded, skipping."
        return 0
    fi
    local zip="data_archive/tred/TRED.sql.zip"
    if [ ! -f "$zip" ]; then
        fail "$zip not found - obtain a copy of the TRED database backup and place it at $zip, then re-run."
    fi
    log "Loading TRED database from $zip"
    local tmp
    tmp=$(mktemp -d)
    unzip -o -q "$zip" -d "$tmp"
    mysql -u"$DB_USER" -p"$DB_PWD" -e 'CREATE DATABASE IF NOT EXISTS TRED'
    mysql -u"$DB_USER" -p"$DB_PWD" TRED < "$tmp"/*.sql
    rm -rf "$tmp"
    log "TRED database loaded."
}

cmd_encode() {
    local dir="$DATA_SET_DIR/encode"
    mkdir -p "$dir"
    if [ -f "$dir/tf-targets.txt" ]; then
        log "$dir/tf-targets.txt already present, skipping."
        return 0
    fi
    if [ -f data_archive/encode/tf-targets.txt ]; then
        cp data_archive/encode/tf-targets.txt "$dir/tf-targets.txt"
        log "Staged the already-archived ENCODE file into $dir (frozen data, per the manual - not re-fetched live)."
        return 0
    fi
    log "data_archive/encode/tf-targets.txt not found - falling back to encode_data_fetch.sh, but"
    log "note its source (archive.gersteinlab.org) is documented as no longer accessible."
    (cd data_archive/encode && ./encode_data_fetch.sh)
    cp data_archive/encode/tf-targets.txt "$dir/tf-targets.txt"
}

cmd_gene_exp_static() {
    mkdir -p "$DATA_SET_DIR/microarray/Pavlidis" "$DATA_SET_DIR/microarray/PrietoCarlos"
    cp -n data_archive/LeeGeneExp/GeneExpWith3FromPavlidis.txt "$DATA_SET_DIR/microarray/Pavlidis/" \
        2>/dev/null || true
    cp -n data_archive/PrietoGeneExp/union60.txt "$DATA_SET_DIR/microarray/PrietoCarlos/" \
        2>/dev/null || true
    log "Static gene-expression files staged (Lee/Prieto - frozen data, never re-downloaded)"
}

# Prints configured-vs-latest-available version for the sources known to drift.
# Always exits 0 - bumping build.params is a human decision, not automatic, per
# the note that some releases (e.g. Pfam without pfamA_interactions.txt.gz) are
# not actually usable even though they're "newer".
cmd_check_latest() {
    log "Pfam: configured=$PFAM_RELEASE"
    local pfam_latest
    pfam_latest=$(curl -fsSL https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/ 2>/dev/null \
        | grep -oE '[0-9]+\.[0-9]+' | sort -V | tail -1 || true)
    if [ -n "${pfam_latest:-}" ]; then
        log "Pfam: latest listed=$pfam_latest"
        if curl -fsIL "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/${pfam_latest}/database_files/pfamA_interactions.txt.gz" \
            >/dev/null 2>&1; then
            log "Pfam: $pfam_latest ships pfamA_interactions.txt.gz - safe to bump PFAM_RELEASE"
        else
            log "Pfam: $pfam_latest does NOT ship pfamA_interactions.txt.gz - keep PFAM_RELEASE=$PFAM_RELEASE"
        fi
    else
        log "Pfam: could not determine latest release (network or listing format changed)"
    fi

    log "BioGRID: configured=$BIOGRID_VERSION"
    local biogrid_latest
    biogrid_latest=$(curl -fsSL https://downloads.thebiogrid.org/BioGRID/Release-Archive/ 2>/dev/null \
        | grep -oE 'BIOGRID-[0-9]+\.[0-9]+\.[0-9]+' | sort -V | tail -1 || true)
    log "BioGRID: latest listed=${biogrid_latest:-unknown}"
    log "BioGRID: verify the new version's tab2/IDENTIFIERS files still match the naming cmd_biogrid"
    log "         expects before bumping (BioGRID has changed major format versions before)."

    log "STRING: configured=$STRING_VERSION"
    # string-db.org's download page doesn't reliably list "the latest" - it can show
    # one version as the default download while banner-advertising a newer one - so
    # this only sanity-checks that the *configured* version is still actually served,
    # rather than guessing at a "latest".
    if curl -fsIL "https://stringdb-downloads.org/download/protein.links.full.${STRING_VERSION}/9606.protein.links.full.${STRING_VERSION}.txt.gz" \
        >/dev/null 2>&1; then
        log "STRING: configured version $STRING_VERSION still resolves."
    else
        log "STRING: configured version $STRING_VERSION did NOT resolve - check https://string-db.org/cgi/download"
    fi
    log "STRING: check https://string-db.org/cgi/download manually for a newer version - the site"
    log "        doesn't expose a clean machine-readable 'latest' listing."

    log "Panther: configured=$PANTHER_VERSION"
    log "Panther: recent releases are BioPAX-only (see manual notes) - verify PantherToReactomeConverterTest"
    log "         still parses the format before bumping PANTHER_VERSION, don't rely on this check alone."
}

# These datasets are frozen/archival by design (the manual explicitly says they
# don't get updated, or the upstream source is known dead) - `all` skips them by
# default so a routine yearly run doesn't re-fetch/re-fail on them. Run any of
# them by name directly, or pass --include-static to `all` to include them.
STATIC_STEPS=(nci-pid tred encode gene-exp-static bioplex)
YEARLY_STEPS=(uniprot pir reactome-dump pfam panther string biogrid go ensembl)

cmd_all() {
    local include_static=0
    if [ "${1:-}" = "--include-static" ]; then
        include_static=1
    fi
    local steps=("${YEARLY_STEPS[@]}")
    if [ "$include_static" -eq 1 ]; then
        steps+=("${STATIC_STEPS[@]}")
    fi
    local results=()
    for step in "${steps[@]}"; do
        if "$0" "$step"; then
            results+=("OK   $step")
        else
            results+=("FAIL $step")
        fi
    done
    echo
    echo "===== fetch_datasets summary ====="
    printf '%s\n' "${results[@]}"
    if [ "$include_static" -eq 0 ]; then
        echo
        echo "Skipped (frozen/static - run once ever, or per-name when actually needed):"
        printf '  %s\n' "${STATIC_STEPS[@]}"
        echo "Re-run with 'all --include-static' to include them, or run any by name directly."
    fi
}

usage() {
    cat <<EOF
Usage: $0 <subcommand>
  Yearly:  uniprot pir reactome-dump pfam panther string biogrid go ensembl
  Static:  nci-pid tred encode gene-exp-static bioplex
           (frozen/archival - not run by 'all' unless --include-static is given)
  check-latest
  all [--include-static]
EOF
}

cmd="${1:-}"
case "$cmd" in
    uniprot) cmd_uniprot ;;
    pir) cmd_pir ;;
    reactome-dump) cmd_reactome_dump ;;
    pfam) cmd_pfam ;;
    panther) cmd_panther ;;
    string) cmd_string ;;
    biogrid) cmd_biogrid ;;
    bioplex) cmd_bioplex ;;
    go) cmd_go ;;
    ensembl) cmd_ensembl ;;
    nci-pid) cmd_nci_pid ;;
    tred) cmd_tred ;;
    encode) cmd_encode ;;
    gene-exp-static) cmd_gene_exp_static ;;
    check-latest) cmd_check_latest ;;
    all) cmd_all "${2:-}" ;;
    *) usage; exit 1 ;;
esac
