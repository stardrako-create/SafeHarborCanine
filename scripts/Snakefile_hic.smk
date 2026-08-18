import os

configfile: "config_hic.yaml"
shell.executable("/bin/bash")
SRA_BIN = "/mnt/d/Biblioteca/Ferramentas Bioinformática/SRA_Toolkit/sratoolkit.3.4.1-win64/bin"
shell.prefix(f'set -euo pipefail; export PATH="{SRA_BIN}:$PATH"; ')

SAMPLES = config["samples"]
SRA_DIR = config["paths"]["sra_dir"]
FASTQ_DIR = config["paths"]["fastq_dir"]
REF_FASTA = config["paths"]["ref_fasta"]
REF_DIR = config["paths"]["ref_dir"]
WORK = config["paths"]["work_dir"]
CHROM_SIZES = os.path.join(REF_DIR, "chrom.sizes")

THREADS = config["params"]["threads_per_lib"]
MIN_MAPQ = config["params"]["min_mapq"]

PER_LIB = os.path.join(WORK, "per_lib")
BWA_INDEX_DONE = REF_FASTA + ".bwt"


rule all:
    input:
        expand(os.path.join(PER_LIB, "{lib}", "pairs", "{lib}.valid.pairs.gz"), lib=SAMPLES),
        expand(os.path.join(PER_LIB, "{lib}", "pairs", "{lib}.dedup.stats.txt"), lib=SAMPLES),


rule sra_download_and_convert:
    """No pre-downloaded .sra for this dataset (unlike ATAC/RRBS) - these 3
    Mischka Hi-C libraries (~35-45GB each) are fetched live via prefetch,
    one at a time (wip=1 caps concurrency so we never have two ~40GB
    downloads or fastq expansions in flight simultaneously)."""
    output:
        r1=temp(os.path.join(FASTQ_DIR, "{lib}", "{lib}_1.fastq")),
        r2=temp(os.path.join(FASTQ_DIR, "{lib}", "{lib}_2.fastq")),
    threads: THREADS
    resources:
        wip=1,
    shell:
        """
        FINALDIR="$(dirname "{output.r1}")"
        mkdir -p "$FINALDIR"
        WORKDIR="{WORK}/sra_convert/{wildcards.lib}"
        rm -rf "$WORKDIR"
        mkdir -p "$WORKDIR/sra" "$WORKDIR/out" "$WORKDIR/tmp"
        prefetch "{wildcards.lib}" -O "$WORKDIR/sra" --max-size u
        fasterq-dump "$WORKDIR/sra/{wildcards.lib}/{wildcards.lib}.sra" --split-files -e {threads} -O "$WORKDIR/out" -t "$WORKDIR/tmp"
        mv "$WORKDIR/out/{wildcards.lib}_1.fastq" "{output.r1}"
        mv "$WORKDIR/out/{wildcards.lib}_2.fastq" "{output.r2}"
        rm -rf "$WORKDIR"
        """


rule align_parse_sort:
    """bwa mem -5SP (chimeric-aware, standard for Hi-C per 4DN/distiller-nf)
    piped straight into pairtools parse + sort, so the ~40-80GB unfiltered
    alignment never touches disk as a BAM/SAM - only the far more compact
    sorted .pairs.gz is written."""
    input:
        r1=rules.sra_download_and_convert.output.r1,
        r2=rules.sra_download_and_convert.output.r2,
        idx=BWA_INDEX_DONE,
    output:
        sorted_pairs=temp(os.path.join(PER_LIB, "{lib}", "pairs", "{lib}.sorted.pairs.gz")),
    threads: THREADS
    resources:
        wip=1,
    shell:
        """
        OUTDIR="$(dirname "{output.sorted_pairs}")"
        TMPDIR="$OUTDIR/tmp_{wildcards.lib}"
        mkdir -p "$OUTDIR" "$TMPDIR"
        bwa mem -5SP -T0 -t {threads} "{REF_FASTA}" "{input.r1}" "{input.r2}" \
            | pairtools parse --chroms-path "{CHROM_SIZES}" --add-columns mapq \
                --nproc-in {threads} --nproc-out {threads} \
            | pairtools sort --tmpdir "$TMPDIR" --nproc {threads} -o "{output.sorted_pairs}"
        rm -rf "$TMPDIR"
        """


rule pairtools_dedup:
    input:
        sorted_pairs=rules.align_parse_sort.output.sorted_pairs,
    output:
        deduped=temp(os.path.join(PER_LIB, "{lib}", "pairs", "{lib}.dedup.pairs.gz")),
        stats=os.path.join(PER_LIB, "{lib}", "pairs", "{lib}.dedup.stats.txt"),
    threads: THREADS
    shell:
        """
        pairtools dedup --nproc-in {threads} --nproc-out {threads} \
            --output-stats "{output.stats}" -o "{output.deduped}" "{input.sorted_pairs}"
        """


rule pairtools_select:
    """Keep only unique-unique pairs above the MAPQ floor - the no-frag
    substitute for exact restriction-fragment filtering, since the Dovetail
    kit's enzyme wasn't disclosed in the source paper."""
    input:
        deduped=rules.pairtools_dedup.output.deduped,
    output:
        valid=os.path.join(PER_LIB, "{lib}", "pairs", "{lib}.valid.pairs.gz"),
    threads: THREADS
    params:
        min_mapq=MIN_MAPQ,
    shell:
        """
        pairtools select \
            '(pair_type == "UU") and (mapq1 >= {params.min_mapq}) and (mapq2 >= {params.min_mapq})' \
            --nproc-in {threads} --nproc-out {threads} \
            -o "{output.valid}" "{input.deduped}"
        """
