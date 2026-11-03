#!/bin/bash
set -eu
out='/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/08_CANDIDATE_CLOSURE_AUDIT/W11_CARRIER_ACCESS_20260929/local_sra_test_20260930'
mkdir -p "$out"
date -u +%FT%TZ > "$out/started.txt"
/usr/bin/fastq-dump-orig --split-files -N 1 -X 100 -A SRR14750332 -O "$out" '/mnt/c/Users/Utilizador/SafeHarborData/SRR14750332/SRR14750332.sra' > "$out/stdout.log" 2> "$out/stderr.log"
date -u +%FT%TZ > "$out/completed.txt"
