import os

configfile: "config_rrbs.yaml"
shell.executable("/bin/bash")
SRA_BIN = "/mnt/d/Biblioteca/Ferramentas Bioinformática/SRA_Toolkit/sratoolkit.3.4.1-win64/bin"
shell.prefix(f'set -euo pipefail; export PATH="{SRA_BIN}:$PATH"; ')

SAMPLES = config["samples"]
SRA_DIR = config["paths"]["sra_dir"]
FASTQ_DIR = config["paths"]["fastq_dir"]
REF_FASTA = config["paths"]["ref_fasta"]
REF_DIR = config["paths"]["ref_dir"]
WORK = config["paths"]["work_dir"]

THREADS = config["params"]["threads_per_sample"]

PER_DOG = os.path.join(WORK, "per_dog")
SCRIPTS_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))
GENOME_PREP_DONE = os.path.join(REF_DIR, "Bisulfite_Genome", "CT_conversion", "BS_CT.1.bt2")


def build_sra_map(sra_dir):
    """SRA download nesting under this dir is inconsistent (some runs sit
    at sra/{ACC}.sra, others at sra/{ACC}/{ACC}.sra, at least one 3 levels
    deep) - walk the tree and match by exact basename instead of assuming a
    fixed depth."""
    mapping = {}
    for root, _dirs, files in os.walk(sra_dir):
        for fname in files:
            if fname.endswith(".sra"):
                mapping[fname[:-4]] = os.path.join(root, fname)
    return mapping


SRA_MAP = build_sra_map(SRA_DIR)


rule all:
    input:
        expand(os.path.join(PER_DOG, "{sample}", "qc", "{sample}.rrbs_qc.tsv"), sample=SAMPLES),
        expand(os.path.join(PER_DOG, "{sample}", "meth", "{sample}_pe.bismark.cov.gz"), sample=SAMPLES),


rule sra_to_fastq:
    input:
        sra=lambda wc: SRA_MAP[wc.sample],
    output:
        r1=temp(os.path.join(FASTQ_DIR, "{sample}", "{sample}_1.fastq")),
        r2=temp(os.path.join(FASTQ_DIR, "{sample}", "{sample}_2.fastq")),
    threads: THREADS
    resources:
        wip=1,
    shell:
        """
        FINALDIR="$(dirname "{output.r1}")"
        mkdir -p "$FINALDIR"
        WORKDIR="{WORK}/sra_convert/{wildcards.sample}"
        rm -rf "$WORKDIR"
        mkdir -p "$WORKDIR/sra" "$WORKDIR/out" "$WORKDIR/tmp"
        ln "{input.sra}" "$WORKDIR/sra/{wildcards.sample}.sra" 2>/dev/null \
            || cp "{input.sra}" "$WORKDIR/sra/{wildcards.sample}.sra"
        fasterq-dump "$WORKDIR/sra/{wildcards.sample}.sra" --split-files -e {threads} -O "$WORKDIR/out" -t "$WORKDIR/tmp"
        mv "$WORKDIR/out/{wildcards.sample}_1.fastq" "{output.r1}"
        mv "$WORKDIR/out/{wildcards.sample}_2.fastq" "{output.r2}"
        rm -rf "$WORKDIR"
        """


rule trim:
    input:
        r1=rules.sra_to_fastq.output.r1,
        r2=rules.sra_to_fastq.output.r2,
    output:
        r1=temp(os.path.join(PER_DOG, "{sample}", "trim", "{sample}_val_1.fq")),
        r2=temp(os.path.join(PER_DOG, "{sample}", "trim", "{sample}_val_2.fq")),
    threads: THREADS
    resources:
        wip=1,
    shell:
        """
        OUTDIR="$(dirname "{output.r1}")"
        mkdir -p "$OUTDIR"
        trim_galore --rrbs --paired --cores {threads} -o "$OUTDIR" \
            --basename "{wildcards.sample}" "{input.r1}" "{input.r2}"
        """


rule bismark_align:
    input:
        r1=rules.trim.output.r1,
        r2=rules.trim.output.r2,
        genome_prep=GENOME_PREP_DONE,
    output:
        bam=temp(os.path.join(PER_DOG, "{sample}", "align", "{sample}_pe.bam")),
        report=os.path.join(PER_DOG, "{sample}", "align", "{sample}_PE_report.txt"),
    threads: THREADS * 2
    resources:
        wip=1,
    params:
        bt_threads=THREADS,
    shell:
        """
        OUTDIR="$(dirname "{output.bam}")"
        TMPDIR="$OUTDIR/tmp_{wildcards.sample}"
        mkdir -p "$OUTDIR" "$TMPDIR"
        bismark --genome "{REF_DIR}" -1 "{input.r1}" -2 "{input.r2}" \
            -p {params.bt_threads} -o "$OUTDIR" -B "{wildcards.sample}" --temp_dir "$TMPDIR"
        rm -rf "$TMPDIR"
        """


rule methylation_extract:
    input:
        bam=rules.bismark_align.output.bam,
    output:
        cov=os.path.join(PER_DOG, "{sample}", "meth", "{sample}_pe.bismark.cov.gz"),
        splitting_report=os.path.join(PER_DOG, "{sample}", "meth", "{sample}_pe_splitting_report.txt"),
    threads: THREADS
    shell:
        """
        OUTDIR="$(dirname "{output.cov}")"
        mkdir -p "$OUTDIR"
        bismark_methylation_extractor -p --comprehensive --bedGraph --cytosine_report \
            --genome_folder "{REF_DIR}" -o "$OUTDIR" "{input.bam}"
        """


rule qc_metrics_rrbs:
    input:
        pe_report=rules.bismark_align.output.report,
        splitting_report=rules.methylation_extract.output.splitting_report,
        cov=rules.methylation_extract.output.cov,
    output:
        qc=os.path.join(PER_DOG, "{sample}", "qc", "{sample}.rrbs_qc.tsv"),
    shell:
        """
        mkdir -p "$(dirname "{output.qc}")"
        python "{SCRIPTS_DIR}/qc_metrics_rrbs.py" \
            --sample "{wildcards.sample}" \
            --pe-report "{input.pe_report}" \
            --splitting-report "{input.splitting_report}" \
            --cov "{input.cov}" \
            --out "{output.qc}"
        """
