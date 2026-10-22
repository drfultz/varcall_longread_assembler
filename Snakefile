localrules: read_splitter, cons_calls, merge_cons

configfile: "config.yaml"

wildcard_constraints:
    contig=r"[^/]+",
    seg_num=r"\d+"

##### Main Assembly Pipeline
#####

#Python script finds matching anchors in feature calls then uses this to
#bin the DNA sequence segments into the appropriate multi-fasta for that repeat unit
checkpoint read_splitter:
    input:
        config["Feature-calls_file"],
        config["Oriented_fastq_file"],
        "loci_anchors/{contig}_anchors.tsv"
    output:
        "split_read_info/{contig}_read_info.tsv",
        directory("segment_reads/{contig}"),
        "metadata/{contig}_metadata.tsv"
    params:
        Segment_q_cutoff=config["Segment_q_cutoff"],
        Max_depth=config["Max_depth"]
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/split_reads_into_units_sm.py"

#kalign2 makes an MSA of each multi-fasta of segment reads from a given unit
rule kalign2:
    input:
        "segment_reads/{contig}/seg_{seg_num}_reads.fa"
    output:
        "MSA_align/{contig}_seg_{seg_num}.aln"
    params:
        ppn="1",
        mem_per_cpu="10gb",
        walltime="24:00:00"
    conda:
        "envs/kalign2.yaml"
    shell: 
        "kalign -gpo 60 -gpe 9 -tgpe 0 -bonus 0 {input} > {output[0]}"

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
        """cons -name {wildcards.contig}_{wildcards.seg_num}_cons -sequence {input[0]} -filter | awk 'NR==1 {{print}} NR>1 {{gsub("[nN]",""); seq=seq $0}} END {{print seq}}' > {output[0]}"""

#Checkpoint function to reevaluate DAG for converging from multiple segments
#back into a single consensus
def merge_cons_input(wildcards):
    seg_dir = checkpoints.read_splitter.get(contig=wildcards.contig).output[1] 
    seg_nums = glob_wildcards(os.path.join(seg_dir, r"seg_{seg_num,\d+}_reads.fa")).seg_num
    return expand("consensus/{contig}_seg_{seg_num}_cons.fa", contig=wildcards.contig, seg_num=sorted(seg_nums, key=int))

# pastes the consensus sequence together in order
rule merge_cons:
    input:
        cons=merge_cons_input
    output:
        "contig_sequence/{contig}_contig.fa"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/merge_cons_sm.py"

