rule filter_mit:
    input:
        expand("{out}/{sample}/combined_with_pipline/coverage_by_mit_histogram.png",out=config["output_dir"],sample=config["sequencing"]),
        expand("{out}/{sample}/combined_with_pipline/mit_reads.fastq",out=config["output_dir"],sample=config["sequencing"]),
        expand("{out}/{sample}/combined_with_pipline/non_mit_reads.fastq",out=config["output_dir"],sample=config["sequencing"])
        
rule map_sbs:
    input:
        ref= config["Mit"],
        sbs="{out}/{sample}/combined_with_pipline/toxo.fastq",
        
    output:
        "{out}/{sample}/combined_with_pipline/sbs_mapped.paf"
    shell:
        "minimap2 -cx map-ont {input.ref} {input.sbs} > {output}"

rule count_sbs_per_read:
    input:
        "{out}/{sample}/combined_with_pipline/sbs_mapped.paf"
    output:
        csv="{out}/{sample}/combined_with_pipline/coverage_by_mit.csv",
        txt="{out}/{sample}/combined_with_pipline/read_names_of_mit_origin.txt",
    shell:
        """
        python ../scripts/calculate_mit_coverage.py {input} {output.csv} {output.txt}
        """
rule histogram:
     input:
          "{out}/{sample}/combined_with_pipline/coverage_by_mit.csv"
     output:
           "{out}/{sample}/combined_with_pipline/coverage_by_mit_histogram.png"
     shell:
           """      
           python ../scripts/histogram.py {input} {output}
           """  
rule filter_fastq:
    input:
        fastq="{out}/{sample}/combined_with_pipline/toxo.fastq",
        ids="{out}/{sample}/combined_with_pipline/read_names_of_mit_origin.txt"
    output:
        "{out}/{sample}/combined_with_pipline/mit_reads.fastq"
    shell:
        "seqkit grep -f {input.ids} {input.fastq} -o {output}"
 
rule filter_non_mitochondrial_reads:
    input:
        fastq="{out}/{sample}/combined_with_pipline/toxo.fastq",
        ids="{out}/{sample}/combined_with_pipline/read_names_of_mit_origin.txt"
    output:
       "{out}/{sample}/combined_with_pipline/non_mit_reads.fastq"
    shell:
        "seqkit grep -v -f {input.ids} {input.fastq} -o {output}"        
                 
