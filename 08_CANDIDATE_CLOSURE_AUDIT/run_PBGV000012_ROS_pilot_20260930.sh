#!/bin/bash
set -euo pipefail
base='/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino'
out="$base/08_CANDIDATE_CLOSURE_AUDIT/SV_PBGV000012/ROS_pilot_20260930"
mkdir -p "$out"
exec 9>"$out/worker.lock"
flock -n 9 || exit 73
trap 'rc=$?; echo "$rc" > "$out/exit_code.txt"; date -u +%FT%TZ > "$out/finished.txt"' EXIT
[ ! -e "$out/pilot.bam" ] || exit 74
date -u +%FT%TZ > "$out/started.txt"
echo extracting > "$out/stage.txt"
/usr/bin/fastq-dump-orig --split-files -N 1 -X 100000 -A SRR14750332 -O "$out" '/mnt/c/Users/Utilizador/SafeHarborData/SRR14750332/SRR14750332.sra' > "$out/extract.stdout" 2> "$out/extract.stderr"
echo aligning > "$out/stage.txt"
bin='/home/stardrako/miniforge3/envs/atac/bin'
"$bin/bwa" mem -t 4 -R '@RG\tID:SRR14750332\tSM:PBGV000012\tPL:ILLUMINA' "$base/01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa" "$out/SRR14750332_1.fastq" "$out/SRR14750332_2.fastq" 2> "$out/bwa.stderr" | "$bin/samtools" sort -@ 2 -m 256M -o "$out/pilot.bam" - 2> "$out/sort.stderr"
"$bin/samtools" quickcheck -v "$out/pilot.bam"
"$bin/samtools" index "$out/pilot.bam"
"$bin/samtools" flagstat "$out/pilot.bam" > "$out/flagstat.txt"
sha256sum "$out/pilot.bam" > "$out/pilot.sha256"
echo completed > "$out/stage.txt"