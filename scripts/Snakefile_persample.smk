import os

configfile: "config.yaml"
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

SAMPLES = config["samples"]
FASTQ_DIR = config["paths"]["fastq_dir"]
REF_FASTA = config["paths"]["ref_fasta"]
REF_GFF = config["paths"]["ref_gff"]
WORK = config["paths"]["work_dir"]

THREADS = config["params"]["threads_per_sample"]
INDEX_THREADS = config["params"]["index_threads"]
GSIZE = config["params"]["effective_genome_size"]
TSS_FLANK = config["params"]["tss_flank"]
TSS_BG = config["params"]["tss_bg_offset"]

REF_DIR = os.path.join(WORK, "reference")
BT2_INDEX = os.path.join(REF_DIR, "ROS_Cfam_1.0")
CHROM_SIZES = os.path.join(REF_DIR, "chrom.sizes")
MITO_LIST = os.path.join(REF_DIR, "mito_contigs.txt")
KEEP_CONTIGS = os.path.join(REF_DIR, "keep_contigs.txt")
TSS_BED = os.path.join(REF_DIR, "tss.bed")

PER_DOG = os.path.join(WORK, "per_dog")
SCRIPTS_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))


rule all:
    input:
        expand(os.path.join(PER_DOG, "{sample}", "qc", "{sample}.qc.tsv"), sample=SAMPLES),
        expand(os.path.join(PER_DOG, "{sample}", "bigwig", "{sample}.cpm.bw"), sample=SAMPLES),
        expand(os.path.join(PER_DOG, "{sample}", "bigwig", "{sample}.raw.bw"), sample=SAMPLES),
        expand(os.path.join(PER_DOG, "{sample}", "peaks", "{sample}_peaks.narrowPeak"), sample=SAMPLES),


rule reference_prep:
    input:
        fasta=REF_FASTA,
        gff=REF_GFF,
    output:
        idx_done=os.path.join(REF_DIR, "bowtie2_index.done"),
        chrom_sizes=CHROM_SIZES,
        mito=MITO_LIST,
        keep=KEEP_CONTIGS,
        tss=TSS_BED,
    threads: INDEX_THREADS
    shell:
        """
        mkdir -p "{REF_DIR}"
        bowtie2-build --threads {threads} "{input.fasta}" "{BT2_INDEX}"
        touch "{output.idx_done}"

        samtools faidx "{input.fasta}"
        cut -f1,2 "{input.fasta}.fai" > "{output.chrom_sizes}"

        grep "^>" "{input.fasta}" | grep -i "mitochondrion" | sed -E 's/^>([^ ]+).*/\\1/' > "{output.mito}" || true

        cut -f1 "{output.chrom_sizes}" | grep -v -x -F -f "{output.mito}" > "{output.keep}"

        python "{SCRIPTS_DIR}/extract_tss.py" --gff "{input.gff}" --out "{output.tss}"
        """


rule trim:
    input:
        r1=os.path.join(FASTQ_DIR, "{sample}", "{sample}_1.fastq"),
        r2=os.path.join(FASTQ_DIR, "{sample}", "{sample}_2.fastq"),
    output:
        r1=temp(os.path.join(PER_DOG, "{sample}", "trim", "{sample}_1.trim.fastq.gz")),
        r2=temp(os.path.join(PER_DOG, "{sample}", "trim", "{sample}_2.trim.fastq.gz")),
        json=os.path.join(PER_DOG, "{sample}", "trim", "{sample}.fastp.json"),
        html=os.path.join(PER_DOG, "{sample}", "trim", "{sample}.fastp.html"),
    threads: THREADS
    shell:
        """
        mkdir -p "$(dirname "{output.r1}")"
        fastp -i "{input.r1}" -I "{input.r2}" -o "{output.r1}" -O "{output.r2}" \
              -j "{output.json}" -h "{output.html}" -w {threads} -q 20
        """


rule align:
    input:
        r1=rules.trim.output.r1,
        r2=rules.trim.output.r2,
        idx_done=os.path.join(REF_DIR, "bowtie2_index.done"),
    output:
        bam=temp(os.path.join(PER_DOG, "{sample}", "align", "{sample}.sorted.bam")),
    threads: THREADS
    shell:
        """
        mkdir -p "$(dirname "{output.bam}")"
        bowtie2 --very-sensitive -X 2000 -p {threads} -x "{BT2_INDEX}" \
            -1 "{input.r1}" -2 "{input.r2}" \
          | samtools sort -@ {threads} -o "{output.bam}" -
        samtools index "{output.bam}"
        """


rule dedup_filter:
    input:
        bam=rules.align.output.bam,
        keep=KEEP_CONTIGS,
    output:
        bam=os.path.join(PER_DOG, "{sample}", "filtered", "{sample}.filt.bam"),
        bai=os.path.join(PER_DOG, "{sample}", "filtered", "{sample}.filt.bam.bai"),
        flagstat=os.path.join(PER_DOG, "{sample}", "filtered", "{sample}.flagstat.txt"),
    threads: THREADS
    shell:
        """
        OUTDIR="$(dirname "{output.bam}")"
        mkdir -p "$OUTDIR"
        TMP_NS="$OUTDIR/{wildcards.sample}.namesort.bam"
        TMP_FM="$OUTDIR/{wildcards.sample}.fixmate.bam"
        TMP_CS="$OUTDIR/{wildcards.sample}.coordsort.bam"
        TMP_MD="$OUTDIR/{wildcards.sample}.markdup.bam"

        samtools sort -n -@ {threads} -o "$TMP_NS" "{input.bam}"
        samtools fixmate -m -@ {threads} "$TMP_NS" "$TMP_FM"
        samtools sort -@ {threads} -o "$TMP_CS" "$TMP_FM"
        samtools markdup -@ {threads} "$TMP_CS" "$TMP_MD"
        samtools index -@ {threads} "$TMP_MD"

        mapfile -t KEEP_ARR < "{input.keep}"
        samtools view -b -q 30 -f 2 -F 1804 -@ {threads} -o "{output.bam}" "$TMP_MD" "${{KEEP_ARR[@]}}"
        samtools index "{output.bam}"
        samtools flagstat "{output.bam}" > "{output.flagstat}"

        rm -f "$TMP_NS" "$TMP_FM" "$TMP_CS" "$TMP_MD"
        """


rule macs2_peaks:
    input:
        bam=rules.dedup_filter.output.bam,
    output:
        narrowpeak=os.path.join(PER_DOG, "{sample}", "peaks", "{sample}_peaks.narrowPeak"),
    params:
        outdir=lambda wc: os.path.join(PER_DOG, wc.sample, "peaks"),
        gsize=GSIZE,
    shell:
        """
        mkdir -p "{params.outdir}"
        macs3 callpeak -t "{input.bam}" -f BAMPE -g {params.gsize} \
            --keep-dup all -q 0.05 -n "{wildcards.sample}" --outdir "{params.outdir}"
        """


rule bigwig:
    input:
        bam=rules.dedup_filter.output.bam,
        bai=rules.dedup_filter.output.bai,
    output:
        bw=os.path.join(PER_DOG, "{sample}", "bigwig", "{sample}.cpm.bw"),
    threads: THREADS
    shell:
        """
        mkdir -p "$(dirname "{output.bw}")"
        bamCoverage -b "{input.bam}" -o "{output.bw}" --normalizeUsing CPM -bs 25 -p {threads}
        """


rule raw_bigwig:
    input:
        bam=rules.dedup_filter.output.bam,
        bai=rules.dedup_filter.output.bai,
    output:
        bw=os.path.join(PER_DOG, "{sample}", "bigwig", "{sample}.raw.bw"),
    threads: THREADS
    shell:
        """
        mkdir -p "$(dirname "{output.bw}")"
        bamCoverage -b "{input.bam}" -o "{output.bw}" -bs 25 -p {threads}
        """


rule qc_metrics:
    input:
        flagstat=rules.dedup_filter.output.flagstat,
        bam=rules.dedup_filter.output.bam,
        peaks=rules.macs2_peaks.output.narrowpeak,
        tss=TSS_BED,
        chrom_sizes=CHROM_SIZES,
    output:
        qc=os.path.join(PER_DOG, "{sample}", "qc", "{sample}.qc.tsv"),
    params:
        tss_flank=TSS_FLANK,
        tss_bg=TSS_BG,
    shell:
        """
        mkdir -p "$(dirname "{output.qc}")"
        python "{SCRIPTS_DIR}/qc_metrics.py" \
            --sample "{wildcards.sample}" \
            --bam "{input.bam}" \
            --peaks "{input.peaks}" \
            --tss "{input.tss}" \
            --chrom-sizes "{input.chrom_sizes}" \
            --tss-flank {params.tss_flank} \
            --tss-bg {params.tss_bg} \
            --out "{output.qc}"
        """
