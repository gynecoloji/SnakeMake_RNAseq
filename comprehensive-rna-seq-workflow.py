import os
import glob

# Collect all sample names from fastq files
SAMPLES, = glob_wildcards("data/{sample}_R1_001.fastq.gz")

rule all:
    input:
        # FastQC and trimming outputs
        expand("results/fastqc/raw/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/raw/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("results/trimmed/{sample}_R1.trimmed.fastq.gz", sample=SAMPLES),
        expand("results/trimmed/{sample}_R2.trimmed.fastq.gz", sample=SAMPLES),
        
        # Alignment outputs
        expand("results/hisat2/{sample}.sam", sample=SAMPLES),
        expand("results/samtools/{sample}.sorted.filtered.bam", sample=SAMPLES),
        expand("results/samtools/{sample}_summary.txt", sample=SAMPLES),
        
        # Salmon quantification
        expand("results/quants/{sample}_quant/quant.sf", sample=SAMPLES),
        expand("results/quants_decoy/{sample}_quant/quant.sf", sample=SAMPLES),
        
        # FeatureCounts quantification
        "results/featurecounts/featureCount.txt",
        
        # RNA-seq QC outputs
        expand("results/picard/{sample}_final_insert_size.txt", sample=SAMPLES),
        expand("results/qualimap_bamqc/{sample}/{sample}.pdf", sample=SAMPLES),
        expand("results/qualimap_rnaseq/{sample}", sample=SAMPLES),
        expand("results/rseqc/{sample}_RD_summary.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}_GC_content.GC.xls", sample=SAMPLES),
        expand("results/rseqc/{sample}.sorted.filtered.tin.xls", sample=SAMPLES),
        
        # Final report
        "results/multiqc_report.html"

# Rule: FastQC on raw reads
rule fastqc_raw:
    input:
        r1="data/{sample}_R1_001.fastq.gz",
        r2="data/{sample}_R2_001.fastq.gz"
    output:
        html_r1="results/fastqc/raw/{sample}_R1_001_fastqc.html",
        html_r2="results/fastqc/raw/{sample}_R2_001_fastqc.html"
    log:
        "logs/fastqc/{sample}.log"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p results/fastqc/raw logs/fastqc
        fastqc -t {threads} -o results/fastqc/raw {input.r1} {input.r2} &> {log}
        """

# Rule: Trim reads with fastp
rule fastp_trim:
    input:
        r1="data/{sample}_R1_001.fastq.gz",
        r2="data/{sample}_R2_001.fastq.gz"
    output:
        r1_trimmed="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2_trimmed="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        html="results/trimmed/{sample}_fastp.html",
        json="results/trimmed/{sample}_fastp.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p results/trimmed logs/fastp
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_trimmed} -O {output.r2_trimmed} \
              --html {output.html} --json {output.json} \
              --thread {threads} \
              --detect_adapter_for_pe -g -p \
              &> {log}
        """

# Rule: FastQC on trimmed reads
rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        html_r1="results/fastqc/trimmed/{sample}_R1.trimmed_fastqc.html",
        html_r2="results/fastqc/trimmed/{sample}_R2.trimmed_fastqc.html"
    log:
        "logs/fastqc/{sample}_trimmed.log"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p results/fastqc/trimmed logs/fastqc
        fastqc -t {threads} -o results/fastqc/trimmed {input.r1} {input.r2} &> {log}
        """

# Rule: Align reads with HISAT2
rule hisat2_align:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        sam="results/hisat2/{sample}.sam",
        summary="results/hisat2/{sample}.sam.summary"
    log:
        "logs/hisat2/{sample}.log"
    params:
        index="ref/ENSEMBL/genome"
    threads: 20
    conda:
        "envs/align.yaml"
    shell:
        """
        mkdir -p results/hisat2 logs/hisat2
        hisat2 -x {params.index} -1 {input.r1} -2 {input.r2} \
               -S {output.sam} \
               --summary-file {output.summary} \
               -p {threads} \
               -q --phred33 -X 3000 -I 0 --no-discordant --no-mixed \
               &> {log}
        """

# Rule: Filter, sort, index BAM and generate flagstat
rule samtools_sort_filter_index:
    input:
        "results/hisat2/{sample}.sam"
    output:
        bam="results/samtools/{sample}.sorted.filtered.bam",
        bai="results/samtools/{sample}.sorted.filtered.bam.bai",
        flagstat="results/samtools/{sample}_summary.txt"
    log:
        "logs/samtools/{sample}.log"
    threads: 20
    conda:
        "envs/align.yaml"
    shell:
        """
        mkdir -p results/samtools logs/samtools
        # Filter for properly paired and uniquely mapped reads
        samtools view -@ {threads} -f 0x2 -F 0x100 -hS {input} | grep "NH:i:1\\|^@" > results/samtools/tmp_{wildcards.sample}.sam
        # Convert to BAM, sort, index and get stats
        samtools view -@ {threads} -bhS results/samtools/tmp_{wildcards.sample}.sam | \
        samtools sort -@ {threads} -O bam -o {output.bam}
        samtools index -@ {threads} {output.bam} {output.bai}
        samtools flagstat {output.bam} > {output.flagstat}
        rm results/samtools/tmp_{wildcards.sample}.sam
        """

# Rule: Salmon quantification (standard index)
rule salmon_quant:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        "results/quants/{sample}_quant/quant.sf"
    log:
        "logs/salmon/{sample}_salmon.log"
    threads: 20
    params:
        index="ref/Salmon_index_Grch38"
    conda:
        "envs/salmon.yaml"
    shell:
        """
        mkdir -p results/quants logs/salmon
        salmon quant \
            -i {params.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            --validateMappings \
            -o results/quants/{wildcards.sample}_quant > {log} 2>&1
        """

# Rule: Salmon quantification (decoy-aware index)
rule salmon_decoy_quant:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        "results/quants_decoy/{sample}_quant/quant.sf"
    log:
        "logs/salmon/{sample}_salmon_decoy.log"
    threads: 20
    params:
        index="ref/Salmon_index_decoy_Grch38"
    conda:
        "envs/salmon.yaml"
    shell:
        """
        mkdir -p results/quants_decoy logs/salmon
        salmon quant \
            -i {params.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            --validateMappings \
            -o results/quants_decoy/{wildcards.sample}_quant > {log} 2>&1
        """

# Rule: Sort BAM by name for paired-end analysis
rule samtools_sort_by_name:
    input:
        "results/samtools/{sample}.sorted.filtered.bam"
    output:
        "results/samtools_byname/{sample}.sorted.byname.bam"
    log:
        "logs/samtools_byname/{sample}.log"
    threads: 8
    conda:
        "envs/align.yaml"
    shell:
        """
        mkdir -p results/samtools_byname logs/samtools_byname
        samtools sort -n -@ {threads} -o {output} {input} &> {log}
        """

# Rule: Quantify with FeatureCounts
rule featurecount:
    input:
        expand("results/samtools/{sample}.sorted.filtered.bam", sample=SAMPLES)
    output:
        "results/featurecounts/featureCount.txt"
    log:
        "logs/featurecounts/featurecount.log"
    params:
        gtf="ref/Homo_sapiens.GRCh38.102.gtf"
    threads: 20
    conda:
        "envs/quant.yaml"
    shell:
        """
        mkdir -p results/featurecounts logs/featurecounts
        featureCounts -T {threads} --countReadPairs -p \
                      -a {params.gtf} -F GTF -t exon -g gene_id -M -s 0 \
                      -o {output} {input} &> {log}
        """

# Rule: Picard CollectInsertSizeMetrics
rule picard_mean_fragment_length:
    input:
        bam="results/samtools/{sample}.sorted.filtered.bam"
    output:
        metrics="results/picard/{sample}_insert_size_metrics.txt",
        final_metrics="results/picard/{sample}_final_insert_size.txt"
    log:
        "logs/picard/{sample}_MeanFragmentLength.log"
    threads: 4
    params:
        picard_jar="ref/picard.jar"
    conda:
        "envs/qualimap.yaml"
    shell:
        """ 
        mkdir -p results/picard logs/picard
        java -jar {params.picard_jar} CollectInsertSizeMetrics \
            I={input.bam} \
            O={output.metrics} \
            H=results/picard/{wildcards.sample}_Histogram.pdf \
            M=0.5 &> {log}
        
        awk 'NR==7,NR==8 {{print $0}}' {output.metrics} > {output.final_metrics}
        """

# Rule: Qualimap BAMQC
rule qualimap_bamqc:
    input:
        bam="results/samtools/{sample}.sorted.filtered.bam"
    output:
        directory("results/qualimap_bamqc/{sample}"),
        pdf="results/qualimap_bamqc/{sample}/{sample}.pdf"
    log:
        "logs/qualimap/{sample}_bamqc.log"
    threads: 8
    params:
        gtf="ref/Homo_sapiens.GRCh38.102.gtf"
    conda:
        "envs/qualimap.yaml"
    shell:
        """
        mkdir -p results/qualimap_bamqc/{wildcards.sample} logs/qualimap
        qualimap bamqc \
            --java-mem-size=20G \
            -bam {input.bam} \
            -c \
            --feature-file {params.gtf} \
            -outdir results/qualimap_bamqc/{wildcards.sample} \
            -os \
            -outfile {wildcards.sample}.pdf \
            -outformat PDF \
            --sequencing-protocol non-strand-specific \
            &> {log}
        """

# Rule: Qualimap RNAseq QC
rule qualimap_rnaseq:
    input:
        bam="results/samtools_byname/{sample}.sorted.byname.bam"
    output:
        directory("results/qualimap_rnaseq/{sample}")
    log:
        "logs/qualimap/{sample}_rnaseq.log"
    threads: 8
    params:
        gtf="ref/Homo_sapiens.GRCh38.102.gtf"
    conda:
        "envs/qualimap.yaml"
    shell:
        """
        mkdir -p results/qualimap_rnaseq/{wildcards.sample} logs/qualimap
        qualimap rnaseq \
            --java-mem-size=20G \
            -a uniquely-mapped-reads \
            -bam {input.bam} \
            -gtf {params.gtf} \
            -outdir results/qualimap_rnaseq/{wildcards.sample} \
            --sequencing-protocol non-strand-specific \
            --paired \
            --sorted \
            &> {log}
        """

# Rule: RSeQC Read distribution
rule rseqc_read_distribution:
    input:
        bam="results/samtools/{sample}.sorted.filtered.bam"
    output:
        "results/rseqc/{sample}_RD_summary.txt"
    log:
        "logs/rseqc/{sample}_read_distribution.log"
    params:
        refbed="ref/ENSEMBL_hg38.bed"
    threads: 4
    conda:
        "envs/RSeQC.yaml"
    shell:
        """
        mkdir -p results/rseqc logs/rseqc
        read_distribution.py -i {input.bam} -r {params.refbed} > {output} 2> {log}
        """

# Rule: RSeQC GC content
rule rseqc_gc_content:
    input:
        bam="results/samtools/{sample}.sorted.filtered.bam"
    output:
        "results/rseqc/{sample}_GC_content.GC.xls"
    log:
        "logs/rseqc/{sample}_gc_content.log"
    threads: 4
    conda:
        "envs/RSeQC.yaml"
    shell:
        """
        mkdir -p results/rseqc logs/rseqc
        read_GC.py -i {input.bam} -o results/rseqc/{wildcards.sample}_GC_content &> {log}
        """

# Rule: RSeQC TIN (Transcript Integrity Number)
rule rseqc_tin:
    input:
        bam="results/samtools/{sample}.sorted.filtered.bam"
    output:
        "results/rseqc/{sample}.sorted.filtered.tin.xls"
    log:
        "logs/rseqc/{sample}_TIN.log"
    params:
        refbed="ref/ENSEMBL_hg38.bed"
    threads: 4
    conda:
        "envs/RSeQC.yaml"
    shell:
        """
        mkdir -p results/rseqc logs/rseqc
        cd results/rseqc && tin.py -i ../../{input.bam} -r ../../{params.refbed} &> ../../{log}
        """

# Rule: Summarize everything with MultiQC
rule multiqc:
    input:
        # FastQC reports
        expand("results/fastqc/raw/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/raw/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/trimmed/{sample}_R1.trimmed_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/trimmed/{sample}_R2.trimmed_fastqc.html", sample=SAMPLES),
        
        # Trimming reports
        expand("results/trimmed/{sample}_fastp.json", sample=SAMPLES),
        
        # Alignment reports
        expand("results/hisat2/{sample}.sam.summary", sample=SAMPLES),
        expand("results/samtools/{sample}_summary.txt", sample=SAMPLES),
        
        # Quantification
        expand("results/quants/{sample}_quant/quant.sf", sample=SAMPLES),
        expand("results/quants_decoy/{sample}_quant/quant.sf", sample=SAMPLES),
        "results/featurecounts/featureCount.txt",
        
        # QC reports
        expand("results/picard/{sample}_insert_size_metrics.txt", sample=SAMPLES),
        expand("results/qualimap_bamqc/{sample}/{sample}.pdf", sample=SAMPLES),
        expand("results/qualimap_rnaseq/{sample}", sample=SAMPLES),
        expand("results/rseqc/{sample}_RD_summary.txt", sample=SAMPLES),
        expand("results/rseqc/{sample}_GC_content.GC.xls", sample=SAMPLES),
        expand("results/rseqc/{sample}.sorted.filtered.tin.xls", sample=SAMPLES)
    output:
        "results/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p logs/multiqc
        multiqc results/ -o results/ -f &> {log}
        """
