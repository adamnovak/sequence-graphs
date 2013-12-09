#!/usr/bin/env bash

# jsonify.sh: take two Avro-format files for a sequence graph (AlleleGroups and
# Adjacencies) and produce a JSON-format file that can be visualized in a
# browser with D3.

# Usage: jsonify.sh test.ag test.adj > test.json

# Die on error
set -e

# We output a JSON object with "allele_groups" and "adjacencies" arrays.
echo "{"
echo "\"allele_groups\": ["

# Make sure to comma-separate the Avro output, without trailing commas. See
# http://stackoverflow.com/a/14518804/402891
avro cat $1 | awk 'NR > 1 { printf(",\n") } { printf($0) }'

echo
echo "], \"adjacencies\": ["

avro cat $2 | awk 'NR > 1 { printf(",\n") } { printf($0) }'

echo
echo "]"
echo "}"
