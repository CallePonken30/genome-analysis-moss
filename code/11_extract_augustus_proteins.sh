#!/bin/bash

INPUT="analyses/03_annotation/augustus_chr3/augustus_chr3_with_prot.gff"
OUTPUT="analyses/03_annotation/augustus_chr3/proteins.fa"

awk '
BEGIN {
    seq=""
    header=""
    inprot=0
}

/^# start gene / {
    gene=$4
    next
}

/^# protein sequence = \[/ {
    header=">" gene
    line=$0
    sub(/^# protein sequence = \[/, "", line)

    if (line ~ /\]$/) {
        sub(/\]$/, "", line)
        print header
        print line
        inprot=0
    } else {
        print header
        print line
        inprot=1
    }
    next
}

inprot {
    line=$0
    sub(/^# /, "", line)

    if (line ~ /\]$/) {
        sub(/\]$/, "", line)
        print line
        inprot=0
    } else {
        print line
    }
}
' "$INPUT" > "$OUTPUT"

echo "Proteins written to $OUTPUT"
