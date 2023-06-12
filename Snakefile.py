SRA,FRR = glob_wildcards("rawreads/{sra}_{frr}.fastq.gz") 

rule rawFastqc:
    input:
        rawread="rawReads/{sra}_{frr}.fastq.gz",
    output:
        zip="rawQC/{sra}_{frr}_fastqc.zip",
        html="rawQC/{sra}_{frr}_fastqc.html",
    threads:
        16
    params:
        path="rawQC/",
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path}
        """
        
        rule get_genome_fa:
    "Downloading Genome sequence, Mus Musculus primary assembly (GRCm39)"
    output:
        fa = "genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
    shell:
        "cd genome"
        " && wget ftp://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
        " && gunzip -k Mus_musculus.GRCm39.dna_sm.primary_assembly.fa "

rule index:
    input:
        fa = rules.get_genome_fa.output.fa
    output:
        dir = ["index."  + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    threads:
        16
    shell:
        " hisat2-build -p {threads} {input.fa} index --quiet"
        
rule hisat_align:
    input:
        fastq1 = rules.fastp.output.read1,
        fastq2 = rules.fastp.output.read2,
        index = rules.index.output.dir
    output:
        bams  = "aligned/{sra}.bam",
        sum   = "logs/{sra}_sum.txt",
        met   = "logs/{sra}_met.txt"
    message:
        "mapping reads to genome to bam files."
    threads: 
        16
    shell:
        "hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x index \
        -1 {input.fastq1} -2 {input.fastq2} | samtools view -Sb -F 4 -o {output.bams}"
