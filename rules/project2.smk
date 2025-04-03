rule project:
    input:
        expand("{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/collective_projected_positions.bed",out=config["output_dir"],sample=config["assemble"])
        
rule assembly_to_fasta2:
     input:
        reads="{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly/40-polishing/filtered_contigs.fasta",
        ref_mmi=config["project"]["fasta"]
     output:
        paf="{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/mapped_reads.paf"
     wildcard_constraints:
        sample="[^/]+"
     threads: 40    
     shell:
        """
        minimap2 -c -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 --secondary=no -t 20 --eqx -Y -O 5,56 -E 4,1 -B 5 {input.ref_mmi} {input.reads} > {output.paf}
        """
rule chain_assembly_to_reference2:
    input:
        paf="{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/mapped_reads.paf"
    output:
        "{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/mapped_reads.chain"
    shell:
        """
        paf2chain -i {input.paf} > {output}
        """

rule liftover2:
    input:
        chain="{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/mapped_reads.chain",
        bed=config["project"]["position"]
    output:
        mapped=temp("{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/projected_positions.bed"),
        unmapped=temp("{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/non_projected_positions.bed"),
        collective="{out}/{sample}/combined_without_pipline/assembly_with_mit_removal/assembly_to_reference/collective_projected_positions.bed"
    shell:
        """
        liftOver -minMatch=0.4 {input.bed} {input.chain} {output.mapped} {output.unmapped}
        python ../scripts/collective_bed.py {input.bed} {output.mapped} {output.unmapped} {output.collective}
        """
      
