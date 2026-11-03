#!/bin/bash
out='/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/08_CANDIDATE_CLOSURE_AUDIT/W11_CARRIER_ACCESS_20260929/local_sra_test_20260930'
exec 9>"$out/validate.lock"
flock -n 9 || exit 73
date -u +%FT%TZ > "$out/validate_started.txt"
vdb-validate -C yes '/mnt/c/Users/Utilizador/SafeHarborData/SRR14750332/SRR14750332.sra' > "$out/validate_stdout.log" 2> "$out/validate_stderr.log"
echo $? > "$out/validate_exit_code.txt"
date -u +%FT%TZ > "$out/validate_finished.txt"
