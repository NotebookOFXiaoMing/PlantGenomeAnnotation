GENOMES, = glob_wildcards(config["ref"])
PROT, = glob_wildcards(config["prot"])

print(GENOMES)
print(PROT)


rule all:
    input:
        expand("04.protein/01.miniprot.gff/{genome}/{prot}.gff3",genome=GENOMES,prot=PROT),
        expand("04.protein/01.miniprot.gff/{genome}/prot.gff3",genome=GENOMES)



rule miniprot:
    input:
        ref = config["ref"],
        prot = config["prot"]
    output:
        "04.protein/01.miniprot.gff/{genome}/{prot}.gff3"
    threads:
        8
    resources:
        mem_mb = 24000
    shell:
        """
        miniprot -t {threads} --gff {input.ref} {input.prot} > {output}
        """


def fixedGenome(wildcards):
    return expand(rules.miniprot.output,genome=wildcards.genome,prot=PROT)

rule cat:
    input:
        fixedGenome
    output:
        "04.protein/01.miniprot.gff/{genome}/prot.gff3"
    threads:
        1
    shell:
        """
        cat {input} > {output}
        """
