GENOMES, = glob_wildcards(config['ref'])

print(GENOMES)

rule all:
    input:
        expand("03.braker2/{genome}/{genome}.braker2/braker.gtf",genome=GENOMES)


rule braker2:
    input:
        fa = config['ref'],
        bam = config['bam']
    output:
        "03.braker2/{genome}/{genome}.braker2/braker.gtf"
    threads:
        32
    resources:
        mem_mb = 64000
    params:
        species = "{genome}",
        output_dir = "03.braker2/{genome}/{genome}.braker2"
    shell:
        """
        time braker.pl --cores {threads} \
        --species={params.species} \
        --softmasking \
        --genome={input.fa} \
        --bam={input.bam} \
        --workingdir={params.output_dir}
        """