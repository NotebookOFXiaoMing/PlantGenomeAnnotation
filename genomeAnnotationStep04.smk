GENOMES, = glob_wildcards(config['ref'])

print(GENOMES)

rule all:
    input:
        expand("05.evm/{genome}/{genome}.pep.fa",genome=GENOMES)


rule augustus:
    input:
        "03.braker2/{genome}/{genome}.braker2/augustus.hints.gtf"
    output:
        "05.evm/{genome}/evm_augustus.gff3"
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl {input} > {output}
        """

rule genemark:
    input:
        "03.braker2/{genome}/{genome}.braker2/GeneMark-ET/genemark.gtf"
    output:
        "05.evm/{genome}/evm_genemark.gff3"
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/misc/GeneMarkHMM_GTF_to_EVM_GFF3.pl {input} > {output}
        """

rule denovo:
    input:
        augustus = rules.augustus.output,
        genemark = rules.genemark.output
    output:
        "05.evm/{genome}/evm_denovo.gff3"
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        cat {input.augustus} {input.genemark} > {output}
        """

rule prot:
    input:
        "04.protein/01.miniprot.gff/{genome}/prot.gff3"
    output:
        "05.evm/{genome}/evm_prot.gff3"
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py {input} > {output}
        """

rule transcript:
    input:
        "02.transcriptome/04.transdecoder/{genome}/transcripts.fa.transdecoder.genome.gff3"
    output:
        "05.evm/{genome}/evm_transcript.gff3"
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        cp {input} {output}
        """




def weightsAbsPath(wildcards):
    import os
    return os.path.abspath("05.evm/" + wildcards.genome + "/weights.txt")


rule partition_EVM_inputs:
    input:
        fa = config["ref"],
        weights = weightsAbsPath,
        denovo = rules.denovo.output,
        prot = rules.prot.output,
        transcript = rules.transcript.output,
    output:
        "05.evm/{genome}/partition_list.out"
    threads:
        32
    resources:
        mem_mb = 80000
    params:
        cd_folder = "05.evm/{genome}/",
        output_folder = "{genome}"
    shell:
        """
        cd {params.cd_folder}
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/partition_EVM_inputs.pl \
        --genome {input.fa} \
        --gene_predictions evm_denovo.gff3 \
        --protein_alignments evm_prot.gff3 \
        --transcript_alignments evm_transcript.gff3 \
        --segmentSize 100000 --overlapSize 10000 \
        --partition_listing partition_list.out \
        --partition_dir {params.output_folder}
        """

rule write_EVM_commands:
    input:
        fa = config["ref"],
        weights = weightsAbsPath,
        denovo = rules.denovo.output,
        prot = rules.prot.output,
        transcript = rules.transcript.output,
        partition_list = rules.partition_EVM_inputs.output
    output:
        "05.evm/{genome}/commands.list"
    threads:
        12
    resources:
        mem_mb = 16000
    params:
        cd_folder = "05.evm/{genome}/",
    shell:
        """
        cd {params.cd_folder}
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/write_EVM_commands.pl \
        --genome {input.fa} \
        --weights {input.weights} \
        --gene_predictions evm_denovo.gff3 \
        --protein_alignments evm_prot.gff3 \
        --transcript_alignments evm_transcript.gff3 \
        --output_file_name evm.out \
        --partitions partition_list.out > commands.list
        """

rule execute_EVM_commands:
    input:
        rules.write_EVM_commands.output
    output:
        "05.evm/{genome}/logs01.txt"
    threads:
        32
    resources:
        mem_mb = 80000
    params:
        "05.evm/{genome}"
    shell:
        """
        cd {params}
        parallel --jobs {threads} < commands.list > logs01.txt
        """

rule recombine_EVM_partial_outputs:
    input:
        rules.execute_EVM_commands.output
    output:
        "05.evm/{genome}/logs02.txt"
    threads:
        12
    resources:
        mem_mb = 24000
    params:
        "05.evm/{genome}"
    shell:
        """
        cd {params}
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/recombine_EVM_partial_outputs.pl \
        --partitions partition_list.out \
        --output_file_name evm.out > logs02.txt
        """

rule convert_EVM_outputs_to_GFF3:
    input:
        log = rules.recombine_EVM_partial_outputs.output,
        fa = config["ref"]
    output:
        "05.evm/{genome}/logs03.txt"
    threads:
        12
    resources:
        mem_mb = 24000
    params:
        "05.evm/{genome}"
    shell:
        """
        cd {params}
        /home/myan/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
        --partitions partition_list.out \
        --output evm.out \
        --genome {input.fa} > logs03.txt
        """

rule find_cat:
    input:
        rules.convert_EVM_outputs_to_GFF3.output
    output:
        "05.evm/{genome}/{genome}_final_EVM.gff3"
    threads:
        12
    resources:
        mem_mb = 24000
    params:
        cd_folder = "05.evm/{genome}",
        gff = "{genome}_final_EVM.gff3"
    shell:
        """
        cd {params.cd_folder}
        find . -regex ".*evm.out.gff3" -exec cat {{}} \\; > {params.gff}
        """

rule generate_key_value:
    input:
        config["ref"]
    output:
        "05.evm/{genome}/key_value.txt"
    threads:
        1
    resources:
        mem_mb = 2000
    params:
        "{genome}"
    run:
        from Bio import SeqIO
        fw = open(output[0],'w')
        i = 0
        for rec in SeqIO.parse(input[0],'fasta'):
            i += 1
            fw.write("%s\t%s\n"%(rec.id,params[0]+str(i).zfill(3)))

        fw.close()

rule rename_gff3:
    input:
        gff = rules.find_cat.output,
        key_value = rules.generate_key_value.output
    output:
        "05.evm/{genome}/{genome}.rename.gff3"
    threads:
        1
    resources:
        mem_mb = 4000
    params:
        folder = "05.evm/{genome}",
        prefix = "{genome}",
        gff = "{genome}_final_EVM.gff3"
    shell:
        """
        cd {params.folder}
        python /data/myan/raw_data/pome/ys.Genome/02.genome.annotation/rename_EVM_output_gff.py -g {params.gff} -c key_value.txt -p {params.prefix}
        """

rule gffread:
    input:
        gff = rules.rename_gff3.output,
        fa = config['ref']
    output:
        cds = "05.evm/{genome}/{genome}.cds.fa",
        pep = "05.evm/{genome}/{genome}.pep.fa"
    threads:
        4
    resources:
        mem_mb = 4000
    shell:
        """
        gffread {input.gff} -g {input.fa} -x {output.cds} -y {output.pep}
        """
