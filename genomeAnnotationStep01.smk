GENOMES, = glob_wildcards(config["ref"]["fa"])
RNASEQ, = glob_wildcards(config["rnaseq"]["r1"])


print(GENOMES)
print(RNASEQ)

rule all:
    input:
        expand("02.transcriptome/02.merged.bam/{genome}/merged.bam.bai",genome=GENOMES),
        expand("02.transcriptome/04.transdecoder/{genome}/transcripts.fa.transdecoder.genome.gff3",genome=GENOMES)


rule fastp:
    input:
        r1 = config['rnaseq']['r1'],
        r2 = config['rnaseq']['r2']
    output:
        r1 = "02.transcriptome/00.clean.reads/{rnaseq}_1.fq",
        r2 = "02.transcriptome/00.clean.reads/{rnaseq}_2.fq",
        html = "02.transcriptome/00.fastp.report/{rnaseq}.html",
        json = "02.transcriptome/00.fastp.report/{rnaseq}.json"
    threads:
        8
    resources:
        mem_mb = 48000
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        -j {output.json} -h {output.html} \
        -w {threads}
        """

rule hisat2_build:
    input:
        fa = config['ref']['fa']
    output:
        "02.transcriptome/00.ref.index/{genome}/ref.1.ht2",
    threads:
        16
    resources:
        mem_mb = 48000
    params:
        "02.transcriptome/00.ref.index/{genome}/ref"
    shell:
        """
        hisat2-build {input} {params}
        """

rule hisat2_align:
    input:
        r1 = rules.fastp.output.r1,
        r2 = rules.fastp.output.r2,
        index = rules.hisat2_build.output
    output:
        "02.transcriptome/01.sam/{genome}/{rnaseq}.sam"
    threads:
        16
    resources:
        mem_mb = 64000
    params:
        "02.transcriptome/00.ref.index/{genome}/ref"
    shell:
        """
        hisat2 -p {threads} --dta -x {params} \
        -1 {input.r1} -2 {input.r2} \
        -S {output}
        """

rule samtools_sort:
    input:
        rules.hisat2_align.output
    output:
        "02.transcriptome/01.sorted.bam/{genome}/{rnaseq}.sorted.bam"
    threads:
        8
    resources:
        mem_mb = 48000
    shell:
        """
        samtools sort -@ {threads} -O bam -o {output} {input}
        """

rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "02.transcriptome/01.sorted.bam/{genome}/{rnaseq}.sorted.bam.bai"
    threads:
        8
    resources:
        mem_mb = 24000
    shell:
        """
        samtools index {input}
        """

def getGenome01(wildcards):
    return expand(rules.samtools_sort.output,genome=wildcards.genome,rnaseq=RNASEQ)

def getGenome02(wildcards):
    return expand(rules.samtools_index.output,genome=wildcards.genome,rnaseq=RNASEQ)

rule samtools_merge:
    input:
        bams = getGenome01,
        bai = getGenome02
    output:
        "02.transcriptome/02.merged.bam/{genome}/merged.bam"
    threads:
        16
    resources:
        mem_mb = 64000
    shell:
        """
        samtools merge -@ {threads} {output} {input.bams}
        """

rule samtools_merge_index:
    input:
        rules.samtools_merge.output
    output:
        "02.transcriptome/02.merged.bam/{genome}/merged.bam.bai"
    threads:
        16
    resources:
        mem_mb = 64000
    shell:
        """
        samtools index {input}
        """

rule stringtie:
    input:
        bam = rules.samtools_merge.output,
        bai = rules.samtools_merge_index.output
    output:
        "02.transcriptome/03.merged.gtf/{genome}/merged.gtf"
    threads:
        24
    resources:
        mem_mb = 48000
    shell:
        """
        stringtie -p {threads} -o {output} {input.bam}
        """


rule transdecoder01:
    input:
        gtf = rules.stringtie.output,
        fa = config['ref']['fa']
    output:
        "02.transcriptome/04.transdecoder/{genome}/transcripts.fa"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        gtf_genome_to_cdna_fasta.pl {input.gtf} {input.fa} > {output}
        """

rule transdecoder02:
    input:
        rules.stringtie.output,
    output:
        "02.transcriptome/04.transdecoder/{genome}/transcripts.gff3"
    threads:
        2
    resources:
        mem_mb = 4000
    shell:
        """
        gtf_to_alignment_gff3.pl {input} > {output}
        """



rule transdecoder03:
    input:
        rules.transdecoder01.output
    output:
        "02.transcriptome/04.transdecoder/{genome}/transcripts.fa.transdecoder_dir/longest_orfs.cds"
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        folder = "02.transcriptome/04.transdecoder/{genome}",
        fa = "transcripts.fa"
    shell:
        """
        cd {params.folder}
        TransDecoder.LongOrfs -t {params.fa}
        """

rule transdecoder04:
    input:
        fa = rules.transdecoder01.output,
        cds = rules.transdecoder03.output
    output:
        "02.transcriptome/04.transdecoder/{genome}/transcripts.fa.transdecoder.gff3"
    threads:
        4
    resources:
        mem_mb = 8000
    params:
        folder = "02.transcriptome/04.transdecoder/{genome}",
        fa = "transcripts.fa"
    shell:
        """
        cd {params.folder}
        TransDecoder.Predict -t {params.fa}
        """

rule transdecoder05:
    input:
        gff3_01 = rules.transdecoder04.output,
        gff3_02 = rules.transdecoder02.output,
        fa = rules.transdecoder01.output
    output:
        "02.transcriptome/04.transdecoder/{genome}/transcripts.fa.transdecoder.genome.gff3"
    threads:
        4
    resources:
        mem_mb = 8000
    shell:
        """
        cdna_alignment_orf_to_genome_orf.pl {input.gff3_01} {input.gff3_02} {input.fa} > {output}
        """