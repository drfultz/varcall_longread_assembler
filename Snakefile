localrules: insert_tags, fq_to_fa, tagged_read_splitter, cons_calls, merge_cons, remove_tags

configfile: "config.yaml"

rule insert_tags:
    input:
        config["Feature-calls_file"],
        config["Oriented_fastq_file"],
        "loci_patterns/{contig}_patterns.tsv",
        config["Tag_set"]
    output:
        "tagged_reads/{contig}_tagreads.fq",
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/insert_tags_multi_sm.py"

rule fq_to_fa:
    input:
        "tagged_reads/{contig}_tagreads.fq"
    output:
        "tagged_reads/{contig}_tagreads.fa"
    params:
        ppn="1",
        vmem="8gb",
        walltime="1:00:00"
    conda:
        "envs/bbmap2.yaml"
    shell:
        "reformat.sh qin=33 in={input[0]} out={output[0]}"

#python script searches reads for sequence in between 'sequence tags',
#then bins the sub-sequences into the appropriate multi-fasta for that tag set
checkpoint tagged_read_splitter:
    input:
        config["Tag_set"],
        "tagged_reads/{contig}_tagreads.fa"
    output:
        directory("segment_reads/{contig}")
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/tagged_read_splitter_sm.py"

#kalign2 makes an MSA of each multi-fasta of reads from a given segment
rule kalign2:
    input:
        "segment_reads/{contig}/seg_{seg_num}_reads.fa"
    output:
        "MSA_align/{contig}_seg_{seg_num}.aln"
    params:
        ppn="1",
        vmem="10gb",
        walltime="24:00:00"
    conda:
        "envs/kalign2.yaml"
    shell: 
        "kalign -gpo 60 -gpe 9 -tgpe 0 -bonus 0 {input} > {output[0]} \n"

#EMBOSS cons makes a consensus of each MSA, followed by removal of 'N' nucleotides
#Default cons settings followed by removal of N is most accurate and
#is robust to low coverage
rule cons_calls:
    input:
        "MSA_align/{contig}_seg_{seg_num}.aln"
    output:
        "consensus/{contig}_seg_{seg_num}_cons.fa"
    conda:
        "envs/cons.yaml"
    shell:
        """cons -name {wildcards.contig}_{wildcards.seg_num}_cons -sequence {input[0]} -filter | awk 'NR<2 {{print}} NR>1 {{gsub("[nN]",""); print}}' > {output[0]}"""
        #"""cons -name {wildcards.contig}_{wildcards.seg_num}_cons -sequence {input[0]} -filter > {output[0]}"""
        #"cons -name {wildcards.contig}_{wildcards.seg_num}_cons -identity 3 -plurality 0.5 -sequence {input[0]} -filter | tr -d 'N' > {output[0]}"  #!!! using tr here will lose 'N's in fa header.

#Checkpoint function to reevaluate DAG for converging from multiple segments
#back into a single consensus
def merge_cons_input(wildcards):
    checkpoint_output = checkpoints.tagged_read_splitter.get(contig=wildcards.contig).output[0]
    return expand("consensus/{contig}_seg_{seg_num}_cons.fa", contig=wildcards.contig, seg_num=glob_wildcards(os.path.join(checkpoint_output, "seg_{seg_num}_reads.fa")).seg_num)

#use echo to write a new fasta header
#then append a single line made from the merged input file sequence lines and
#after removing tandem 'sequence tag' duplications via python script
rule merge_cons:
    input:
        cons=merge_cons_input,
        tag_fa=config["Tag_set"]
    output:
        "merge_cons/{contig}_merge.fa"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/merge_cons_sm.py"
        
#polishing of the draft consensus via minimap2 and racon
#draft and reads still contain tags to ensure alignment at proper locations
rule minimap2_racon:
    input:
        "merge_cons/{contig}_merge.fa",
        "tagged_reads/{contig}_tagreads.fa"
    output:
        "mappings_polish/{contig}_to_raw.sam",
        "polished/{contig}_polished1X.fa"
    params:
        ppn="4",
        vmem="12gb",
        walltime="20:00:00"
    conda:
        "envs/minimap2_racon.yaml"
    shell:
        "minimap2 -a -L --secondary=no -t {params.ppn} -x map-ont {input[0]} {input[1]} | samtools view -h -F 2304 - > {output[0]} \n"
        "racon -t {params.ppn} {input[1]} {output[0]} {input[0]} > {output[1]} \n"

#second round of polishing. Final output still contains tags.
#Visualization of resulting alignments is important to verify proper polishing
rule minimap2_racon2:
    input:
        "polished/{contig}_polished1X.fa",
        "tagged_reads/{contig}_tagreads.fa"
    output:
        "mappings_polish/{contig}_to_1X.sam",
        "polished/{contig}_polished2X.fa"
    params:
        ppn="4",
        vmem="12gb",
        walltime="20:00:00"
    conda:
        "envs/minimap2_racon.yaml"
    shell:
        "minimap2 -a -L --secondary=no -t {params.ppn} -x map-ont {input[0]} {input[1]} | samtools view -h -F 2304 - > {output[0]} \n"
        "racon -t {params.ppn} {input[1]} {output[0]} {input[0]} > {output[1]} \n"

rule remove_tags:
    input:
        config["Tag_set"],
        "polished/{contig}_polished2X.fa"
    output:
        "final_contigs/{contig}_contig.fa"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/remove_tags_sm.py"





#not currently used
rule sizeselect_reads:
    input:
        "",
        ""
    output:
        "reads_sized/{contig}_sizedreads.fq"
    params:
        ppn="1",
        vmem="8gb",
        walltime="2:00:00"
    shell:
        "MEAN=$(awk 'NR % 2 == 0 {{sum+=length}} END {{ print sum / NR * 2 }}' {input[0]} ) \n"
        "let LOWER=MEAN*9/10 \n"
        "let UPPER=MEAN*11/10 \n"
        "awk 'BEGIN {{OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= '$LOWER' && length(seq) <= '$UPPER') {{print header, seq, qheader, qseq}}}}' < {input[1]} > {output}"
