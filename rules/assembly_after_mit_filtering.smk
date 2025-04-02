rule assemble_after_mit_removal: 
     input:
           expand("{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/reads_to_assembly/mapped_reads_sorted.bam",out=config["output_dir"],sample=config["assemble"]),           
rule assembly_2:
    input:
        filter_fastq="{out}/{sample}/combined_with_pipline/non_mit_reads.fastq"
    output:
        assembly=directory("{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/assembly/"),
        out_file="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/assembly/40-polishing/filtered_contigs.fasta",
    wildcard_constraints:
        sample="[^/]+", 
    params:
        genome_size="70m",
    threads: 10
    shell:
        """
        flye --nano-corr {input.filter_fastq} --threads {threads} -o {output.assembly}
        """
rule minimap2_reads_to_assembly_2:
    input:
        ref="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/assembly/40-polishing/filtered_contigs.fasta",
        reads="{out}/{sample}/combined_with_pipline/non_mit_reads.fastq"
    output:
        sam="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/reads_to_assembly/mapped_reads.sam",
        paf="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/reads_to_assembly/mapped_reads.paf"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        """
        minimap2 -ax map-ont {input.ref} {input.reads} > {output.sam}
        minimap2 -cx map-ont {input.ref} {input.reads} > {output.paf}
        """
rule samtools_reads_to_assembly_2:
    input:
        sam="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/reads_to_assembly/mapped_reads.sam"
    output:
        bam="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/reads_to_assembly/mapped_reads.bam",
        sorted="{out}/{sample}/combined_with_pipline/assembly_with_mit_removal/reads_to_assembly/mapped_reads_sorted.bam",
    wildcard_constraints:
        sample="[^/]+"  
    shell:
        """
        samtools view -bS {input.sam} -o {output.bam}
        samtools sort {output.bam} -o {output.sorted}
        samtools index {output.sorted}
        """   
