# Snakemake workflow for the Sunflower Variant Calling Pipeline.
#
# Covers steps 03-06 (trim/align -> HaplotypeCaller -> GenotypeGVCFs -> imputation).
# Steps 00-02 remain as external scripts; run one of:
#   Path A (NCBI):        00_prep_manifest.pbs -> 01_download_srr.pbs -> 02_fastq_convert_concat.pbs
#   Path B (local FASTQs): 00_prep_local_fastq.sh
# Both produce GENLIST.txt and a FASTQ directory as input to this workflow.
#
# Usage:
#   snakemake --profile profiles/pbs    # PBS/Torque cluster
#   snakemake --profile profiles/slurm  # SLURM cluster
#   snakemake --cores N                 # local execution (testing only)
#
# Before running, edit config/config.yaml to set your paths and parameters.
#
# Shell template notation:
#   {wildcards.x}, {input.x}, {params.x}, {config[key]} -- Snakemake substitutions
#   ${{BASH_VAR}} renders as ${BASH_VAR} in the generated shell (bash variable reference)

import os
import glob as _glob

configfile: "config/config.yaml"

# ── Derived reference paths ───────────────────────────────────────────────────
# Mirrors the variables derived from GENOME_DIR in site_config.sh.
_GENOME_DIR = config["genome_dir"]
REF         = os.path.join(_GENOME_DIR, "Ha412HOv2_w_CPMT.fa")
CHROMS_FILE = os.path.join(_GENOME_DIR, "17_chromosomes_Ha412HOv2.list")
REPEAT_BED  = os.path.join(_GENOME_DIR, "helianthusannuus.filteredRepeats.17_chromosomes.bed")

# ── Derived software paths ────────────────────────────────────────────────────
# Mirrors the variables derived from SOFTWARE_DIR in site_config.sh.
_SW      = config["software_dir"]
FASTP    = os.path.join(_SW, "fastp")
BWA_MEM2 = os.path.join(_SW, "bwa-mem2", "bwa-mem2")
MOSDEPTH = os.path.join(_SW, "mosdepth")

# ── Lists read at parse time ──────────────────────────────────────────────────
COHORT  = config["cohort"]
SAMPLES = [s.strip() for s in open(config["genlist"])  if s.strip()]
CHROMS  = [c.strip() for c in open(CHROMS_FILE)        if c.strip()]
POPS    = config["populations"]

# ── FASTQ input helpers ───────────────────────────────────────────────────────
# Replicates the glob from 03_trim_align_markdup.pbs: ${sample}_*1.f*q.gz
# Fails loudly at DAG-construction time if a sample's files are missing or
# ambiguous, rather than silently mid-run.

def _find_fastq(sample, suffix_glob):
    pattern = os.path.join(config["fq_dir"], f"{sample}_{suffix_glob}")
    hits = _glob.glob(pattern)
    if len(hits) == 0:
        raise ValueError(
            f"No FASTQ found for sample '{sample}'.\n"
            f"  Pattern tried: {pattern}\n"
            f"  Check fq_dir in config.yaml and that the sample name matches exactly."
        )
    if len(hits) > 1:
        raise ValueError(
            f"Multiple FASTQs matched for sample '{sample}': {hits}\n"
            f"  Each sample must have exactly one R1 and one R2 in fq_dir."
        )
    return hits[0]

def get_r1(wildcards):
    return _find_fastq(wildcards.sample, "*1.f*q.gz")

def get_r2(wildcards):
    return _find_fastq(wildcards.sample, "*2.f*q.gz")

# ── Pre-flight ────────────────────────────────────────────────────────────────
onstart:
    for d in [
        "logs/trim_align_markdup",
        "logs/haplotype_caller",
        "logs/genotype_gvcfs",
        "logs/imputation",
        "logs/cluster",
    ]:
        os.makedirs(d, exist_ok=True)

# ── Target ────────────────────────────────────────────────────────────────────
localrules: all

rule all:
    input:
        expand(
            os.path.join(config["imputed_dir"], "{cohort}.{pop}.imputed.biallelic_snps.vcf.gz"),
            cohort=COHORT,
            pop=POPS,
        )

# ─────────────────────────────────────────────────────────────────────────────
# Step 03 -- Trim, align, mark duplicates
# One job per sample. Mirrors 03_trim_align_markdup.pbs.
# ─────────────────────────────────────────────────────────────────────────────
rule trim_align_markdup:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        bam = os.path.join(config["bam_dir"], "{sample}.bam"),
        bai = os.path.join(config["bam_dir"], "{sample}.bam.bai"),
    params:
        fastp    = FASTP,
        bwa      = BWA_MEM2,
        mosdepth = MOSDEPTH,
        ref      = REF,
    log:
        "logs/trim_align_markdup/{sample}.log"
    threads: 32
    resources:
        cpus      = 32,
        mem_mb    = 256000,
        walltime  = "168:00:00",
        queue     = config["pbs_queue_large"],
        partition = config["slurm_partition_large"],
    shell:
        """
        exec > {log} 2>&1
        [ -n "{config[modules][samtools]}" ] && module load {config[modules][samtools]}

        mkdir -p {config[bam_dir]}/fastp_out    {config[bam_dir]}/markdup_out \
                 {config[bam_dir]}/mosdepth_out {config[bam_dir]}/flagstat_out

        {params.fastp} \
            -i {input.r1} \
            -I {input.r2} \
            -h {config[bam_dir]}/fastp_out/{wildcards.sample}.html \
            -j {config[bam_dir]}/fastp_out/{wildcards.sample}.json \
            -w 16 \
            --stdout \
        | {params.bwa} \
            mem -p -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" -t 16 {params.ref} - \
        | samtools collate -@ 16 -O - \
        | samtools fixmate  -@ 16 -m - - \
        | samtools sort     -m 7G -@ 16 -T {config[scratch_dir]}/{wildcards.sample} - \
        | samtools markdup  -@ 16 \
            -f {config[bam_dir]}/markdup_out/{wildcards.sample}.stats \
            - {output.bam}

        samtools index {output.bam}
        {params.mosdepth} -t 16 --no-per-base \
            {config[bam_dir]}/mosdepth_out/{wildcards.sample} {output.bam}
        samtools flagstat -@ 16 {output.bam} \
            > {config[bam_dir]}/flagstat_out/{wildcards.sample}.flagstat
        """

# ─────────────────────────────────────────────────────────────────────────────
# Step 04 -- HaplotypeCaller (GVCF mode)
# One job per (sample, chromosome). The PBS version ran all 17 chromosomes in
# one job via GNU parallel; splitting gives finer retry granularity at the same
# total resource cost. The `parallel` module is no longer needed.
# ─────────────────────────────────────────────────────────────────────────────
rule haplotype_caller:
    input:
        bam = os.path.join(config["bam_dir"], "{sample}.bam"),
        bai = os.path.join(config["bam_dir"], "{sample}.bam.bai"),
    output:
        gvcf = os.path.join(config["gvcf_dir"], "{chrom}", "{sample}.{chrom}.g.vcf.gz"),
        tbi  = os.path.join(config["gvcf_dir"], "{chrom}", "{sample}.{chrom}.g.vcf.gz.tbi"),
    params:
        ref        = REF,
        repeat_bed = REPEAT_BED,
    log:
        "logs/haplotype_caller/{sample}.{chrom}.log"
    threads: 1
    resources:
        cpus      = 1,
        mem_mb    = 4000,
        walltime  = "24:00:00",
        queue     = config["pbs_queue_default"],
        partition = config["slurm_partition_default"],
    shell:
        """
        exec > {log} 2>&1
        [ -n "{config[modules][gatk]}" ] && module load {config[modules][gatk]}

        mkdir -p {config[gvcf_dir]}/{wildcards.chrom}
        gatk HaplotypeCaller \
            -R  {params.ref} \
            -I  {input.bam} \
            --intervals               {wildcards.chrom} \
            --exclude-intervals       {params.repeat_bed} \
            --emit-ref-confidence     GVCF \
            --heterozygosity          0.01 \
            --heterozygosity-stdev    0.1 \
            --native-pair-hmm-threads 1 \
            -O  {output.gvcf}
        """

# ─────────────────────────────────────────────────────────────────────────────
# Step 05 -- GenomicsDBImport + GenotypeGVCFs
# One job per chromosome; all samples genotyped jointly.
# Mirrors 05_genotype_gvcfs.pbs.
# ─────────────────────────────────────────────────────────────────────────────
rule genotype_gvcfs:
    input:
        gvcfs = expand(
            os.path.join(config["gvcf_dir"], "{{chrom}}", "{sample}.{{chrom}}.g.vcf.gz"),
            sample=SAMPLES,
        ),
    output:
        vcf = os.path.join(config["vcf_dir"], f"{COHORT}.{{chrom}}.vcf.gz"),
        tbi = os.path.join(config["vcf_dir"], f"{COHORT}.{{chrom}}.vcf.gz.tbi"),
    params:
        ref = REF,
    log:
        f"logs/genotype_gvcfs/{COHORT}.{{chrom}}.log",
    threads: 16
    resources:
        cpus      = 16,
        mem_mb    = 128000,
        walltime  = "168:00:00",
        queue     = config["pbs_queue_default"],
        partition = config["slurm_partition_default"],
    shell:
        """
        exec > {log} 2>&1
        [ -n "{config[modules][gatk]}" ] && module load {config[modules][gatk]}

        mkdir -p {config[vcf_dir]}

        # Build two-column GATK sample map from the explicit GVCF input list.
        # Using Snakemake's input list avoids the `find` call in the PBS version.
        SAMPLE_MAP={config[scratch_dir]}/{config[cohort]}.{wildcards.chrom}.sample_map
        > "${{SAMPLE_MAP}}"
        for gvcf in {input.gvcfs}; do
            sample=$(basename "$gvcf" .{wildcards.chrom}.g.vcf.gz)
            printf '%s\t%s\n' "$sample" "$gvcf" >> "${{SAMPLE_MAP}}"
        done

        gatk --java-options "-Xmx125g -XX:ParallelGCThreads=16" GenomicsDBImport \
            --genomicsdb-workspace-path {config[scratch_dir]}/{config[cohort]}.{wildcards.chrom}.gdb \
            --sample-name-map           "${{SAMPLE_MAP}}" \
            -L                          {wildcards.chrom} \
            --tmp-dir                   {config[scratch_dir]} \
            --reader-threads            16 \
            --overwrite-existing-genomicsdb-workspace

        gatk --java-options "-Xmx125g -XX:ParallelGCThreads=16" GenotypeGVCFs \
            -R  {params.ref} \
            -V  "gendb://{config[scratch_dir]}/{config[cohort]}.{wildcards.chrom}.gdb" \
            -L  {wildcards.chrom} \
            --heterozygosity        0.01 \
            --heterozygosity-stdev  0.1 \
            -O  {output.vcf}
        """

# ─────────────────────────────────────────────────────────────────────────────
# Step 06 -- Imputation and filtering (per population)
# One job per population. Mirrors 06_imputation.pbs.
# ─────────────────────────────────────────────────────────────────────────────
rule imputation:
    input:
        vcfs = expand(
            os.path.join(config["vcf_dir"], "{cohort}.{chrom}.vcf.gz"),
            cohort=COHORT,
            chrom=CHROMS,
        ),
        pop_list = lambda wc: os.path.join(config["pop_list_dir"], f"{wc.pop}.txt"),
    output:
        vcf = os.path.join(config["imputed_dir"], f"{COHORT}.{{pop}}.imputed.biallelic_snps.vcf.gz"),
        tbi = os.path.join(config["imputed_dir"], f"{COHORT}.{{pop}}.imputed.biallelic_snps.vcf.gz.tbi"),
    params:
        mem_gb    = 120,
        n_threads = 16,
    log:
        f"logs/imputation/{COHORT}.{{pop}}.log",
    threads: 16
    resources:
        cpus      = 16,
        mem_mb    = 128000,
        walltime  = "168:00:00",
        queue     = config["pbs_queue_default"],
        partition = config["slurm_partition_default"],
    shell:
        """
        exec > {log} 2>&1
        [ -n "{config[modules][bcftools]}" ] && module load {config[modules][bcftools]}

        mkdir -p {config[imputed_dir]}

        # Concatenate all chromosome VCFs into a genome-wide VCF in scratch.
        CONCAT_VCF={config[scratch_dir]}/{config[cohort]}.{wildcards.pop}.all_chroms.vcf.gz
        bcftools concat \
            --threads {params.n_threads} \
            -a \
            -Oz -o "${{CONCAT_VCF}}" \
            {input.vcfs}
        bcftools index -t "${{CONCAT_VCF}}"

        # Subset to this population's samples.
        POP_VCF={config[scratch_dir]}/{config[cohort]}.{wildcards.pop}.raw.vcf.gz
        bcftools view -S {input.pop_list} \
            --threads {params.n_threads} \
            -Oz -o "${{POP_VCF}}" "${{CONCAT_VCF}}"
        bcftools index -t "${{POP_VCF}}"

        # BEAGLE 5.5 phasing and imputation.
        IMPUTED_PREFIX={config[imputed_dir]}/{config[cohort]}.{wildcards.pop}.imputed
        java -Xmx{params.mem_gb}g -jar {config[beagle_jar]} \
            gt="${{POP_VCF}}" \
            out="${{IMPUTED_PREFIX}}" \
            nthreads={params.n_threads}

        # Filter to biallelic SNPs.
        bcftools view \
            -m2 -M2 -v snps \
            --threads {params.n_threads} \
            -Oz -o {output.vcf} "${{IMPUTED_PREFIX}}.vcf.gz"
        bcftools index -t {output.vcf}

        # Clean up large intermediates.
        rm -f "${{CONCAT_VCF}}" "${{CONCAT_VCF}}.tbi"
        rm -f "${{POP_VCF}}" "${{POP_VCF}}.tbi"
        rm -f "${{IMPUTED_PREFIX}}.vcf.gz" "${{IMPUTED_PREFIX}}.vcf.gz.tbi" \
              "${{IMPUTED_PREFIX}}.log"
        """
