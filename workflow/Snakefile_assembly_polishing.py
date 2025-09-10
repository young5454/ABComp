import os
import glob
# 081924 : Split snakefile - assembly & polishing

configfile: "config/config.yml"

# Main working directory
workspace = config["workspace"]

# Define GROUP, STRAIN and REF wildcards
GROUP, STRAIN = glob_wildcards(os.path.join(workspace, "0.Assembly/{group}_{strain}/genome"))
# REF, = glob_wildcards(os.path.join(workspace, "0.Assembly/ref/genome/{ref}.fasta"))


# Print Message
print("+------------------------------------------+")
print("🧬 ABComp-Assembly-Polishing Run Starting 🐍")
print("+------------------------------------------+\n")

print("📂 GROUPs and STRAINs Detected :")
group_strain_pairs = sorted(zip(GROUP, STRAIN), key=lambda x: x[0])
for g, s in group_strain_pairs:
    print(f" - GROUP : {g: <15} | STRAIN : {s}")

print("+-------------------------------------------+\n")

# Rule all
rule all:
    input:
        ### FastQC-raw I/O files
        # Raw Illumina paired-end short reads
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/reads/{group}_{strain}_1.fastq.gz"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/reads/{group}_{strain}_2.fastq.gz"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # FastQC htmls
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/raw/{group}_{strain}_1_fastqc.html"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/raw/{group}_{strain}_2_fastqc.html"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # FastQC zips
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/raw/{group}_{strain}_1_fastqc.zip"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/raw/{group}_{strain}_2_fastqc.zip"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),   
        ### Trimmomatic I/O files
        # Trimmed Illumina paired-end short reads : Paired
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_1_paired.fastq.gz"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_2_paired.fastq.gz"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # Trimmed Illumina paired-end short reads : Unpaired
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_1_unpaired.fastq.gz"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_2_unpaired.fastq.gz"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        ### FastQC-trimmed I/O files
        # FastQC htmls
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/trimmed/{group}_{strain}_1_paired_fastqc.html"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/trimmed/{group}_{strain}_2_paired_fastqc.html"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # FastQC zips
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/trimmed/{group}_{strain}_1_paired_fastqc.zip"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/trimmed/{group}_{strain}_2_paired_fastqc.zip"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        ### Polypolish I/O files
        # Long-read draft assembly
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/genome/{group}_{strain}.fasta"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # Raw alignments
        expand(
            os.path.join(workspace, 
            "0.Assembly/{group}_{strain}/genome/{group}_{strain}_aligned_1.sam"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "0.Assembly/{group}_{strain}/genome/{group}_{strain}_aligned_2.sam"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # Filtered alignments
        expand(
            os.path.join(workspace,
            "0.Assembly/{group}_{strain}/genome/{group}_{strain}_filtered_1.sam"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "0.Assembly/{group}_{strain}/genome/{group}_{strain}_filtered_2.sam"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # Polished assembly
        expand(
            os.path.join(workspace,
            "0.Assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        ### Busco directory
        expand(
            os.path.join(workspace, "1.Busco/{group}_{strain}_busco/"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        ### Quast directory
        expand(
            os.path.join(workspace, "1.Quast/{group}_{strain}_quast/"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        ### MultiQC dummy
        os.path.join(workspace, "0.QC/multiqc_raw/.multiqc_done"),
        os.path.join(workspace, "0.QC/multiqc_trimmed/.multiqc_done")


# Rule to quality check raw Illumina paired-end short reads
rule fastqc_raw:
    input:
        reads1=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/reads/{group}_{strain}_1.fastq.gz"),
        reads2=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/reads/{group}_{strain}_2.fastq.gz")
    output:
        html1=os.path.join(workspace,
        "0.QC/{group}_{strain}/raw/{group}_{strain}_1_fastqc.html"),
        html2=os.path.join(workspace,
        "0.QC/{group}_{strain}/raw/{group}_{strain}_2_fastqc.html"),
        zip1=os.path.join(workspace,
        "0.QC/{group}_{strain}/raw/{group}_{strain}_1_fastqc.zip"),
        zip2=os.path.join(workspace,
        "0.QC/{group}_{strain}/raw/{group}_{strain}_2_fastqc.zip"),
    message:
        "Run FastQC to raw Illumina paired-end short reads"
    params:
        threads=config["fastqc_threads"],
        out_path=os.path.join(workspace, "0.QC/{group}_{strain}/raw/"),
        multiqc_path=os.path.join(workspace, "0.QC/multiqc_raw/"),
    conda:
        "env_qc"
    shell:
        """
        fastqc \
            -t {params.threads} \
            -o {params.out_path} \
            -q \
            {input.reads1} {input.reads2}
        
        mkdir -p {params.multiqc_path}
        cp {params.out_path}* {params.multiqc_path}
        """


# UPDATE : Rule to integrate raw FastQC results
rule multiqc_raw:
    input:
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/raw/{group}_{strain}_1_fastqc.html"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "0.QC/{group}_{strain}/raw/{group}_{strain}_2_fastqc.html"),
            zip,
            group=GROUP,
            strain=STRAIN,
        )
    output:
        # Dummy output file to ensure rule execution
        os.path.join(workspace, "0.QC/multiqc_raw/.multiqc_done")
    message:
        "Run MultiQC to integrate raw FastQC reports"
    params:
        multiqc_path=os.path.join(workspace, "0.QC/multiqc_raw/")
    conda:
        "env_qc"
    shell:
        """
        multiqc {params.multiqc_path} \
            --force \
            -n fastqc_raw \
            -o {params.multiqc_path}

        touch {output}
        """


# Rule to trim Illumina paired-end short reads
rule trimmomatic:
    input:
        reads1=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/reads/{group}_{strain}_1.fastq.gz"),
        reads2=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/reads/{group}_{strain}_2.fastq.gz")
    output:
        paired1=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_1_paired.fastq.gz"),
        paired2=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_2_paired.fastq.gz"),
        unpaired1=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_1_unpaired.fastq.gz"),
        unpaired2=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/trimmed/{group}_{strain}_2_unpaired.fastq.gz"),
    message:
        "Run Trimmomatic to trim and remove adapters from raw Illumina paired-end short reads"
    params:
        threads=config["trimmomatic_threads"],
        phred=config["trimmomatic_phred"],
        adapter=os.path.join(workspace, config["trimmomatic_adapter"]),
        seed_mismatch=config["trimmomatic_seed_mismatch"],
        palindrome_clip_thres=config["trimmomatic_palindrome_clip_thres"],
        simple_clip_thres=config["trimmomatic_simple_clip_thres"],
        leading=config["trimmomatic_leading"],
        trailing=config["trimmomatic_trailing"],
        sliding1=config["trimmomatic_sliding1"],
        sliding2=config["trimmomatic_sliding2"],
        minlen=config["trimmomatic_minlen"]
    conda:
        "env_trimmomatic"
    shell:
        """
        trimmomatic PE \
            -quiet \
            -threads {params.threads} {params.phred} \
            {input.reads1} {input.reads2} \
            {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} \
            ILLUMINACLIP:{params.adapter}:{params.seed_mismatch}:{params.palindrome_clip_thres}:{params.simple_clip_thres} \
            LEADING:{params.leading} TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.sliding1}:{params.sliding2} MINLEN:{params.minlen}
        """


# Rule to quality check trimmed Illumina paired-end short reads
rule fastqc_trimmed:
    input:
        reads1=rules.trimmomatic.output.paired1,
        reads2=rules.trimmomatic.output.paired2
    output:
        html1=os.path.join(workspace,
        "0.QC/{group}_{strain}/trimmed/{group}_{strain}_1_paired_fastqc.html"),
        html2=os.path.join(workspace,
        "0.QC/{group}_{strain}/trimmed/{group}_{strain}_2_paired_fastqc.html"),
        zip1=os.path.join(workspace,
        "0.QC/{group}_{strain}/trimmed/{group}_{strain}_1_paired_fastqc.zip"),
        zip2=os.path.join(workspace,
        "0.QC/{group}_{strain}/trimmed/{group}_{strain}_2_paired_fastqc.zip"),
        done=os.path.join(workspace, 
        "0.QC/{group}_{strain}/trimmed/.fastqc_trimmed_done")
    message:
        "Run FastQC to trimmed Illumina paired-end short reads"
    params:
        threads=config["fastqc_threads"],
        out_path=os.path.join(workspace, "0.QC/{group}_{strain}/trimmed/"),
        multiqc_path=os.path.join(workspace, "0.QC/multiqc_trimmed/"),
    conda:
        "env_qc"
    shell:
        """
        fastqc \
            -t {params.threads} \
            -o {params.out_path} \
            -q \
            {input.reads1} {input.reads2}

        mkdir -p {params.multiqc_path}
        cp {params.out_path}* {params.multiqc_path}

        touch {output.done}
        """


# UPDATE : Rule to integrate trimmed FastQC results
rule multiqc_trimmed:
    input:
        expand(
            os.path.join(workspace, 
            "0.QC/{group}_{strain}/trimmed/.fastqc_trimmed_done"),
            zip,
            group=GROUP,
            strain=STRAIN
        )
    output:
        # Dummy output file to ensure rule execution
        os.path.join(workspace, "0.QC/multiqc_trimmed/.multiqc_done")
    message:
        "Run MultiQC to integrate trimmed FastQC reports"
    params:
        multiqc_path=os.path.join(workspace, "0.QC/multiqc_trimmed/")
    conda:
        "env_qc"
    shell:
        """
        multiqc {params.multiqc_path} \
            --force \
            -n fastqc_trimmed \
            -o {params.multiqc_path}

        touch {output}
        """


# Rule to run Polypolish for polishing long-read assemblies with Illumina paired-end short reads
rule polypolish:
    input:
        genome_fasta=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}.fasta"),
        reads1=rules.trimmomatic.output.paired1,
        reads2=rules.trimmomatic.output.paired2
    output:
        aligned1=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}_aligned_1.sam"),
        aligned2=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}_aligned_2.sam"),
        filtered1=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}_filtered_1.sam"),
        filtered2=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}_filtered_2.sam"),
        polished_fasta=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"),
    message:
        "Run Polypolish for polishing long-read assemblies with Illumina paired-end short reads"
    log:
        logfile=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/polypolish.log")
    params:
        threads=config["polypolish_threads"],  # Number of threads to use for multi-threaded processes
        path=os.path.join(workspace, "0.Assembly/{group}_{strain}/genome/"),
    conda:
        "env_polypolish"
    shell:
        """
        bwa index {input.genome_fasta} &>> {log.logfile}
        bwa mem -t {params.threads} -a {input.genome_fasta} {input.reads1} > {output.aligned1} 2>> {log.logfile}
        bwa mem -t {params.threads} -a {input.genome_fasta} {input.reads2} > {output.aligned2} 2>> {log.logfile}

        polypolish_insert_filter.py \
            --in1 {output.aligned1} \
            --in2 {output.aligned2} \
            --out1 {output.filtered1} \
            --out2 {output.filtered2} &>> {log.logfile}

        polypolish {input.genome_fasta} {output.filtered1} {output.filtered2} > {output.polished_fasta} 2>> {log.logfile}
        """


# Rule to run Busco assessments for polished assemblies and plot assessment plots
rule busco:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta=rules.polypolish.output.polished_fasta,
    output:
        out_dir=directory(os.path.join(workspace,
        "1.Busco/{group}_{strain}_busco/")),
    message:
        "Run Busco assessments for polished assemblies and plot assessment plots"
    params:
        lineage_path=os.path.join(workspace,
        "busco_downloads/lineages/", config["busco_lineage"]),
        output="{group}_{strain}_busco",
        out_path=os.path.join(workspace, "1.Busco/"),
        sum_dir=os.path.join(workspace, "1.Busco/summaries/"),
    conda:
        "env_busco"
    shell:
        """
        busco \
            -m genome \
            -i {input.strain_fasta} \
            -o {params.output} \
            --out_path {params.out_path} \
            -l {params.lineage_path}

        mkdir -p {params.sum_dir}
        cd {output.out_dir}
        cp *.txt {params.sum_dir}
        
        generate_plot.py -wd {params.sum_dir}
        """


# Rule to run Quast assessments for polished assemblies and plot assessment plots
rule quast:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta=rules.polypolish.output.polished_fasta,
    output:
        out_dir=directory(os.path.join(workspace, "1.Quast/{group}_{strain}_quast/")),
    message:
        "Run Quast assessments for polished assemblies and plot assessment plots"
    params:
        threads=config["quast_threads"],
    conda:
        "env_quast"
    shell:
        """
        quast \
            -o {output.out_dir} \
            --threads {params.threads} {input.strain_fasta}
        """