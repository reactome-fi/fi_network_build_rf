#!/usr/bin/env bash
# Compiles the project directly with javac/a fixed classpath instead of
# `mvn compile`, and writes the resulting runtime classpath to
# target/runtime_classpath.txt for scripts/run_pipeline.sh to use.
#
# Why: on this machine, `mvn compile` fails outright - one of the project's own
# dependencies (org.reactome.fi:modeling:1.0.3) declares a <repositories> entry
# ("drmaa") pointing at a long-dead UMass host, and Maven's dependency resolver
# refuses to trust the already-present local jars without revalidating against
# it (even with -o/offline). The jars pom.xml declares, plus the transitive
# runtime dependencies discovered by actually running the pipeline, are listed
# below explicitly. If you later fix the Maven setup, this script becomes
# unnecessary - `mvn compile exec:java` (what run_pipeline.sh used originally)
# would work again.
#
# If you add a new dependency to pom.xml, add its jar path here too.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

M2="${HOME}/.m2/repository"

DEPS=(
    "$M2/org/reactome/curator-tool/1.0.0/curator-tool-1.0.0.jar"
    "$M2/org/reactome/pathway-exchange/1.0.0/pathway-exchange-1.0.0.jar"
    "$M2/org/reactome/fi/modeling/1.0.3/modeling-1.0.3.jar"
    "$M2/org/reactome/fi/foundation/1.0.3/foundation-1.0.3.jar"
    "$M2/nz/ac/waikato/cms/weka/weka-stable/3.6.6/weka-stable-3.6.6.jar"
    "$M2/log4j/log4j/1.2.17/log4j-1.2.17.jar"
    "$M2/org/slf4j/slf4j-api/2.0.0/slf4j-api-2.0.0.jar"
    "$M2/org/jdom/jdom/1.1.3/jdom-1.1.3.jar"
    "$M2/jaxen/jaxen/1.2.0/jaxen-1.2.0.jar"
    "$M2/org/hibernate/hibernate-core/3.6.10.Final/hibernate-core-3.6.10.Final.jar"
    "$M2/javassist/javassist/3.12.1.GA/javassist-3.12.1.GA.jar"
    "$M2/c3p0/c3p0/0.9.1.2/c3p0-0.9.1.2.jar"
    "$M2/javax/xml/bind/jaxb-api/2.3.1/jaxb-api-2.3.1.jar"
    "$M2/com/sun/xml/bind/jaxb-impl/2.3.4/jaxb-impl-2.3.4.jar"
    "$M2/junit/junit/4.13.2/junit-4.13.2.jar"
    "$M2/org/protege/protege/4.0.0/protege-4.0.0.jar"
    "$M2/org/protege/protege-owl/4.0.0/protege-owl-4.0.0.jar"
    "$M2/org/protege/jena/4.0.0/jena-4.0.0.jar"
    "$M2/org/protege/rdf-api/2001-01-19/rdf-api-2001-01-19.jar"
    "$M2/org/protege/owlsyntax/4.0.0/owlsyntax-4.0.0.jar"
    "$M2/org/protege/xercesImpl/4.0.0/xercesImpl-4.0.0.jar"
    "$M2/tech/tablesaw/tablesaw-core/0.38.4/tablesaw-core-0.38.4.jar"
    "$M2/tech/tablesaw/tablesaw-jsplot/0.38.4/tablesaw-jsplot-0.38.4.jar"
    "$M2/mysql/mysql-connector-java/5.1.47/mysql-connector-java-5.1.47.jar"
    # Transitive runtime dependencies not declared directly in pom.xml, found by
    # actually running the pipeline and adding jars until NoClassDefFoundError
    # stopped appearing:
    "$M2/commons-logging/commons-logging/1.2/commons-logging-1.2.jar"
    "$M2/com/ibm/icu/icu4j/70.1/icu4j-70.1.jar"
    "$M2/dom4j/dom4j/1.6.1/dom4j-1.6.1.jar"
    "$M2/org/hibernate/hibernate-commons-annotations/3.2.0.Final/hibernate-commons-annotations-3.2.0.Final.jar"
    "$M2/antlr/antlr/2.7.7/antlr-2.7.7.jar"
    "$M2/commons-collections/commons-collections/3.2.2/commons-collections-3.2.2.jar"
    "$M2/org/apache/geronimo/specs/geronimo-jta_1.1_spec/1.1.1/geronimo-jta_1.1_spec-1.1.1.jar"
    "$M2/org/hibernate/javax/persistence/hibernate-jpa-2.0-api/1.0.1.Final/hibernate-jpa-2.0-api-1.0.1.Final.jar"
    "$M2/jonelo/jacksum/1.0.0/jacksum-1.0.0.jar"
    "$M2/colt/colt/1.2.0/colt-1.2.0.jar"
    "$M2/com/google/guava/guava/30.0-jre/guava-30.0-jre.jar"
    "$M2/org/apache/commons/commons-math/2.2/commons-math-2.2.jar"
    "$M2/org/jgrapht/jgrapht-core/0.9.1/jgrapht-core-0.9.1.jar"
    "$M2/it/unimi/dsi/fastutil/8.3.0/fastutil-8.3.0.jar"
    # Rest of tablesaw-core 0.38.4's real dependency list (read directly from its
    # own pom.xml in the local repo, rather than discovering each one by crash):
    "$M2/org/apache/commons/commons-math3/3.6.1/commons-math3-3.6.1.jar"
    "$M2/org/roaringbitmap/RoaringBitmap/0.8.12/RoaringBitmap-0.8.12.jar"
    "$M2/com/univocity/univocity-parsers/2.8.4/univocity-parsers-2.8.4.jar"
    "$M2/io/github/classgraph/classgraph/4.8.60/classgraph-4.8.60.jar"
    # javax.activation was removed from the JDK in Java 9+; jaxb-impl needs it
    # at runtime (JAXBContext -> RuntimeModelBuilder -> RuntimeBuiltinLeafInfoImpl)
    # even though it's not needed at compile time.
    "$M2/com/sun/activation/javax.activation/1.2.0/javax.activation-1.2.0.jar"
)

for jar in "${DEPS[@]}"; do
    [ -f "$jar" ] || { echo "Missing dependency jar: $jar" >&2; exit 1; }
done

DEPS_CP=$(IFS=:; echo "${DEPS[*]}")
mkdir -p target/classes
SOURCES_FILE="target/sources.txt"
find src -name "*.java" > "$SOURCES_FILE"
echo "Compiling $(wc -l < "$SOURCES_FILE" | tr -d ' ') source files..."
javac -proc:none -cp "target/classes:resources:$DEPS_CP" -d target/classes "@$SOURCES_FILE"

FULL_CP="target/classes:resources:$DEPS_CP"
echo "$FULL_CP" > target/runtime_classpath.txt
echo "Compiled OK. Runtime classpath written to target/runtime_classpath.txt"
