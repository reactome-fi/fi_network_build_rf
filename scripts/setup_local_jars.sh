#!/usr/bin/env bash
# One-time (per machine, not per year) setup: installs the protege-family jars
# from lib/ into the local Maven repo, using the exact groupId/artifactId/version
# triples already declared in pom.xml. See README.md.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

./install_jar/installJar.sh lib/protege.jar            org.protege protege      4.0.0
./install_jar/installJar.sh lib/protege-owl.jar        org.protege protege-owl  4.0.0
./install_jar/installJar.sh lib/jena.jar               org.protege jena         4.0.0
./install_jar/installJar.sh lib/rdf-api-2001-01-19.jar org.protege rdf-api      2001-01-19
./install_jar/installJar.sh lib/owlsyntax.jar          org.protege owlsyntax    4.0.0
./install_jar/installJar.sh lib/xercesImpl.jar         org.protege xercesImpl   4.0.0

echo "Local protege-family jars installed into the local Maven repo."
echo "Also run the ant scripts in ant/ once (curator-tool, pathway-exchange, jacksum, this project's own jar)."
