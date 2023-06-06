GENOMES, = glob_wildcards(config["ref"])

print(GENOMES)

rule all:
    input:
        expand("05.evm/{genome}/weights.txt",genome=GENOMES)




rule weights:
    input:
        config['weights']
    output:
        "05.evm/{genome}/weights.txt"
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        cp {input} {output}
        """