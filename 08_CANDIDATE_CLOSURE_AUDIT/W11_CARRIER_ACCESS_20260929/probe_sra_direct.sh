#!/bin/bash
set -u
out='/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/08_CANDIDATE_CLOSURE_AUDIT/W11_CARRIER_ACCESS_20260929/sra_probe_direct'
mkdir -p "$out"
date -u +%FT%TZ > "$out/started.txt"
timeout 55s fastq-dump-orig --split-files -N 1 -X 100 -A SRR14750332 -O "$out" 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14750332/SRR14750332' > "$out/stdout.log" 2> "$out/stderr.log"
rc=$?
printf '%s\n' "$rc" > "$out/exit_code.txt"
date -u +%FT%TZ > "$out/finished.txt"
cat "$out/stdout.log" "$out/stderr.log"
exit "$rc"