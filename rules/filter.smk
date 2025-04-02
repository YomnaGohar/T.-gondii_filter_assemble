ruleorder: combine_fastq > statistics > generate_filtered_fastq
out=config["output_dir"]
import os
def generate_input_dict(config):
    input_dict = {}  # Initialize the dictionary
    for sample_wildcards in config["sequencing"].keys():
        input_dict[sample_wildcards] = {}
        for f in config["sequencing"][sample_wildcards]:
            parts = f.split('/')
            sample_directory = '/'.join(parts[:-1])               
            input_dict[sample_directory] = f
    return input_dict
input=generate_input_dict(config)    
def combined_output(assemblies,name):
    output = []  # Initialize the output list
    for f in assemblies[name]:
        parts = f.split('/')
        output_directory = '/'.join(parts[:-1])
        output_path = os.path.join(output_directory,"repeat/filtered_with_pipeline/toxo.fastq")
        output.append(output_path)
    return output
def extract_sample_output(config):
    output = []  # Initialize the output list
    for sample_wildcards in config["sequencing"].keys():
        for f in config["sequencing"][sample_wildcards]:
            parts = f.split('/')
            sample_directory = '/'.join(parts[:-1]) 
            output_path = os.path.join(sample_directory,"repeat/filtered_with_pipeline/toxo.fastq")
            output_path_stat = os.path.join(sample_directory,"repeat/filtered_with_pipeline/stats.HM.txt")
            output.append(output_path)
            output.append(output_path_stat)   
    return output
rule filter:
    input:
         extract_sample_output(config),
         expand("{out}/{sample}/combined_with_pipline/toxo.fastq",out=config["output_dir"],sample=config["sequencing"])
rule minimap2_mapping:
    input:
        fastq=lambda wildcards: input[wildcards.dir],
        minimizer=config["Minimizer"]
    output:
        mapped="{dir}/repeat/filtered_with_pipeline/mapped_sample_to_all_strains_and_mouse_human.bam"
    threads: 10   
    shell:
        """
        minimap2 -ax map-ont {input.minimizer} {input.fastq} -t 10 | samtools view -bS > {output.mapped}
        """

rule convert_bam_to_sam:
    input:
        "{dir}/repeat/filtered_with_pipeline/mapped_sample_to_all_strains_and_mouse_human.bam"
    output:
        "{dir}/repeat/filtered_with_pipeline/noSupp_mapped_sample_to_all_strains_and_mouse_human.sam"
    threads: 5    
    shell:
        """
        samtools view -h -F0x900 {input} > {output}
        """

rule generate_paf:
    input:
        "{dir}/repeat/filtered_with_pipeline/noSupp_mapped_sample_to_all_strains_and_mouse_human.sam"
    output:
        "{dir}/repeat/filtered_with_pipeline/noSupp_mapped_sample_to_all_strains_and_mouse_human.paf"            
    shell:
        """
        export PATH=/home/yogah100/bin/:$PATH
        k8 /gpfs/project/projects/DiltheyLab/projects/Toxoplasma/software/paftools.js sam2paf {input} > {output}
        """

rule filter_human_mouse:
    input:
        "{dir}/repeat/filtered_with_pipeline/noSupp_mapped_sample_to_all_strains_and_mouse_human.paf"
    output:
        "{dir}/repeat/filtered_with_pipeline/names_of_reads_mapped_to_HumanMouse.txt"
    shell:
        """
        awk '$6 ~ /^(NT_|NW_|NC_|chr)/' {input} | awk '{{print $1}}' > {output}
        """

rule filter_not_human_mouse:
    input:
        "{dir}/repeat/filtered_with_pipeline/noSupp_mapped_sample_to_all_strains_and_mouse_human.paf"
    output:
        "{dir}/repeat/filtered_with_pipeline/names_of_reads_NOT_mapped_to_HumanMouse.txt"
    shell:
        """
        awk '$6 !~ /^(NT_|NW_|NC_|chr)/' {input} | awk '{{print $1}}' > {output}
        """

rule generate_filtered_fastq:
    input:
        fastq=lambda wildcards: input[wildcards.dir],
        hm="{dir}/repeat/filtered_with_pipeline/names_of_reads_mapped_to_HumanMouse.txt",
        toxo="{dir}/repeat/filtered_with_pipeline/names_of_reads_NOT_mapped_to_HumanMouse.txt"
    output:
        hm_fastq="{dir}/repeat/filtered_with_pipeline/HM.fastq",
        filter_fastq="{dir}/repeat/filtered_with_pipeline/toxo.fastq",
    shell:
        """
        module load SeqKit
        seqkit grep -f {input.hm} {input.fastq} -o {output.hm_fastq}
        seqkit grep -f {input.hm} -v {input.fastq} -o {output.filter_fastq}
        """
rule statistics:
    input:
        hm=rules.generate_filtered_fastq.output.hm_fastq,
        toxo=rules.generate_filtered_fastq.output.filter_fastq
    output:
        hm_fastq="{dir}/repeat/filtered_with_pipeline/stats.HM.txt",
        filter_fastq="{dir}/repeat/filtered_with_pipeline/stats.toxo.txt"        
    shell:
        """
        module load SeqKit/0.15.0
        seqkit stats -aa {input.hm} > {output.hm_fastq}
        seqkit stats -aa {input.toxo} > {output.filter_fastq}
        """
rule combine_fastq:
    input:
         filter_fastq=lambda wildcards: combined_output(config["sequencing"],wildcards.sample)
    output:
        combined_fastq=f"{out}/{{sample}}/combined_with_pipline/toxo.fastq"
    shell:
        """
        cat {input.filter_fastq}> {output.combined_fastq}
        """    
