localrules: read_splitter, cons_calls, merge_cons, bam_splitter, faidx, merge_mods

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
#back into a single consensus. Shared by both merge steps
# so the contig and the merged BED use the same units in the same order.
def unit_seg_nums(contig):
    seg_dir = checkpoints.read_splitter.get(contig=contig).output[1]
    nums = glob_wildcards(os.path.join(seg_dir, r"seg_{seg_num,\d+}_reads.fa")).seg_num
    return sorted(nums, key=int)

def merge_cons_input(wc):
    return expand("consensus/{contig}_seg_{seg_num}_cons.fa",
                  contig=wc.contig, seg_num=unit_seg_nums(wc.contig))

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




##### Post-assembly Base Modification
#####

#Segment unaligned reads from a sequencing .bam file
checkpoint bam_splitter:
    input:
        config["Bam_reads_file"],
        "split_read_info/{contig}_read_info.tsv",
        "metadata/{contig}_metadata.tsv"
    output:
        directory("segment_bams/{contig}")
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/segment_readbams_sm.py"

#Align the bam segments to their unit's consensus with dorado align (minimap2)
rule dorado_align:
    input:
        "segment_bams/{contig}/seg_{seg_num}_bamreads.bam",
        "consensus/{contig}_seg_{seg_num}_cons.fa"
    output:
        "dorado_align/{contig}_seg{seg_num}_bamread-to-cons.bam",
        "dorado_align/{contig}_seg{seg_num}_bamread-to-cons.bam.bai"
    params:
        ppn="1",
        mem_per_cpu="8gb",
        walltime="1:00:00"
    conda:
        "envs/dorado.yaml"
    shell:
        "dorado aligner --mm2-preset map-ont {input[1]} {input[0]} | samtools sort -o {output[0]} && samtools index {output[0]}"

# Index the fasta for modkit
rule faidx:
    input:
        "consensus/{contig}_seg_{seg_num}_cons.fa"
    output:
        "consensus/{contig}_seg_{seg_num}_cons.fa.fai"
    conda:
        "envs/dorado.yaml"
    shell:
        "samtools faidx {input}"

#Call modifications for each segment
rule modkit:
    input:
        "dorado_align/{contig}_seg{seg_num}_bamread-to-cons.bam",
        "dorado_align/{contig}_seg{seg_num}_bamread-to-cons.bam.bai",
        "consensus/{contig}_seg_{seg_num}_cons.fa",
        "consensus/{contig}_seg_{seg_num}_cons.fa.fai"
    output:
        "modkit_beds/{contig}_seg{seg_num}_modkit.bed"
    params:
        ppn="1",
        mem_per_cpu="8gb",
        walltime="1:00:00",
        flags=config["Modkit_flags"]
    conda:
        "envs/modkit.yaml"
    shell: # Segments with no reads aligned (from either poor mapping or no input reads at dorado_align) are left empty
        """
        if [ "$(samtools view -c -F 4 {input[0]})" -eq 0 ]; then
            touch {output}
        else
            modkit pileup {input[0]} {output} --ref {input[2]} {params.flags}
        fi
        """

### This block makes sure that each consensus segment has a matching modkit bed before proceeding
def merge_mods_input(wc):
    checkpoints.bam_splitter.get(contig=wc.contig)  # make sure BAM segments exist
    segs = unit_seg_nums(wc.contig)
    return {
        "cons": expand("consensus/{contig}_seg_{seg_num}_cons.fa", contig=wc.contig, seg_num=segs),
        "beds": expand("modkit_beds/{contig}_seg{seg_num}_modkit.bed", contig=wc.contig, seg_num=segs),
    }

#Merge the modkit outputs by sorting them and successively
#renumbering the nt positions
rule merge_mods:
    input:
        unpack(merge_mods_input)
    output:
        "modkit_merged/{contig}_merged.bed"
    script:
        "scripts/merge_mods_sm.py"
