import os
import glob
# 082324 : Split snakefile - downstream

configfile: "config/config.yml"

# Main working directory
workspace = config["workspace"]

# Define GROUP, STRAIN and REF wildcards
GROUP, STRAIN = glob_wildcards(os.path.join(workspace, "0.Assembly/{group}_{strain}/genome"))
REF, = glob_wildcards(os.path.join(workspace, "0.Assembly/ref/genome/{ref}.fasta"))

# Print message
print("+--------------------------------------------+")
print("🧬 ABComp-Comparative-Genomics Run Starting 🐍")
print("+--------------------------------------------+\n")

print("📂 GROUPs and STRAINs Detected :")
group_strain_pairs = sorted(zip(GROUP, STRAIN), key=lambda x: x[0])
for g, s in group_strain_pairs:
    print(f" - GROUP : {g: <15} | STRAIN : {s}")

print("\n📌 REFERENCE Genome:")
print(f" - Ref : {REF}\n")

print("+--------------------------------------------+\n")


# Rule all
rule all:
    input:
        ### Prokka I/O files
        # Reference genome & genbank
        expand(os.path.join(workspace, "0.Assembly/ref/genome/{ref}.fasta"), ref=REF),
        # Prokka outputs - REFERENCE 
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.err"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.ffn"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.fsa"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.gff"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.sqn"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.tsv"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.faa"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.fna"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.gbk"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.log"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.tbl"), ref=REF),
        expand(os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.txt"), ref=REF),

        # Prokka outputs - STRAINS 
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.err"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.ffn"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.fsa"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.gff"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.sqn"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.tsv"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.faa"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.fna"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.gbk"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.log"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.tbl"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "2.Annotation/{group}_{strain}/{group}_{strain}.txt"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),

        ### Roary I/O files
        # Roary outputs - STRAIN-to-REF pairwise
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory.header.embl"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace, 
            "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory.tab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory_binary_genes.fa"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory_graph.dot"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/blast_identity_frequency.Rtab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/clustered_proteins"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/core_accessory.header.embl"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/core_accessory.tab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/core_accessory_graph.dot"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/core_alignment_header.embl"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/core_gene_alignment.aln"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/gene_presence_absence.Rtab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/gene_presence_absence.csv"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_conserved_genes.Rtab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_genes_in_pan_genome.Rtab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_new_genes.Rtab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_unique_genes.Rtab"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/pan_genome_reference.fa"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        expand(
            os.path.join(workspace,
            "3.RoaryRef/{group}_{strain}_ref_pairwise/summary_statistics.txt"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
        # Roary output - Within-Group TMP
        expand(os.path.join(workspace, "4.RoaryGroup.tmp/{group}"), group=GROUP),

        # Roary directory - Within-Group Roary
        expand(os.path.join(workspace, "4.RoaryGroup/{group}/"), group=GROUP),

        ### All Group directory
        expand(os.path.join(workspace, "5.AllGroups/{group}/gene_lists/"), group=GROUP),
        expand(os.path.join(workspace, "5.AllGroups/{group}/faas/"), group=GROUP),
        expand(os.path.join(workspace, "5.AllGroups/{group}/FASTA/"), group=GROUP),

        # Roary plots directory
        expand(os.path.join(workspace, "5.AllGroups/{group}/roary_plots/"), group=GROUP),

        ### EggNog-Mapper COG analysis output
        expand(os.path.join(workspace, "5.AllGroups/{group}/emapper-core/"), group=GROUP),
        expand(os.path.join(workspace, "5.AllGroups/{group}/emapper-shells/"), group=GROUP),
        expand(os.path.join(workspace, "5.AllGroups/{group}/cog_plots/"), group=GROUP),

        ### Abricate directory - STRAIN
        expand(os.path.join(workspace,
        "6.Abricate/{group}_{strain}/"), zip, group=GROUP, strain=STRAIN),

        ### Abricate directory - REF
        expand(os.path.join(workspace, "6.Abricate/ref/{ref}"), ref=REF)


# Rule to run Prokka for annotating REFERENCE genome FASTA(s)
rule prokka_ref:
    input:
        ref_fasta=os.path.join(workspace, "0.Assembly/ref/genome/{ref}.fasta"),
        ref_genbank=os.path.join(workspace, "0.Assembly/ref/{ref}.gb"),
    output:
        ref_err=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.err"),
        ref_ffn=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.ffn"),
        ref_fsa=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.fsa"),
        ref_gff=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.gff"),
        ref_sqn=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.sqn"),
        ref_tsv=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.tsv"),
        ref_faa=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.faa"),
        ref_fna=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.fna"),
        ref_gbk=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.gbk"),
        ref_log=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.log"),
        ref_tbl=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.tbl"),
        ref_txt=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.txt"),
    message:
        "Run Prokka for annotating reference genome FASTA(s)"
    params:
        threads=config["prokka_ref_threads"],
        out_dir=os.path.join(workspace, "2.Annotation/ref/{ref}"),
        prefix="{ref}",
        locustag="{ref}",
        kingdom=config["prokka_ref_kingdom"],
    conda:
        "env_prokka"
    shell:
        """
        prokka \
            --outdir {params.out_dir} \
            --prefix {params.prefix} \
            --locustag {params.locustag} \
            --cpus {params.threads} \
            --kingdom {params.kingdom} \
            --addgenes \
            --force \
            --quiet \
            --proteins {input.ref_genbank} {input.ref_fasta}
        """


# Rule to run Prokka for annotating polished STRAIN FASTAs
rule prokka_strain:
    input:
        # Input FASTAs are polished FASTAs from Polypolish
        strain_fasta=os.path.join(workspace,
        "0.Assembly/{group}_{strain}/genome/{group}_{strain}_polished.fasta"),
        strain_ref_genbank=expand(os.path.join(workspace,
        "0.Assembly/ref/{ref}.gb"), ref=REF),
    output:
        strain_err=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.err"),
        strain_ffn=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.ffn"),
        strain_fsa=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.fsa"),
        strain_gff=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.gff"),
        strain_sqn=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.sqn"),
        strain_tsv=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.tsv"),
        strain_faa=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.faa"),
        strain_fna=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.fna"),
        strain_gbk=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.gbk"),
        strain_log=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.log"),
        strain_tbl=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.tbl"),
        strain_txt=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.txt"),
    message:
        "Run Prokka for annotating polished strain genome FASTAs"
    params:
        threads=config["prokka_strain_threads"],
        out_dir=os.path.join(workspace, "2.Annotation/{group}_{strain}/"),
        prefix="{group}_{strain}",
        locustag="{group}_{strain}",
        kingdom=config["prokka_strain_kingdom"],
    conda:
        "env_prokka"
    shell:
        """
        prokka \
            --outdir {params.out_dir} \
            --prefix {params.prefix} \
            --locustag {params.locustag} \
            --cpus {params.threads} \
            --kingdom {params.kingdom} \
            --addgenes \
            --force \
            --quiet \
            --proteins {input.strain_ref_genbank} {input.strain_fasta}
        """


# Rule to run Roary for STRAIN-to-REF 1:1 pairwise Pangenome analysis
rule roary_strain_ref_pairwise:
    input:
        # Input GFFs are annotated GFFs from Prokka
        ref_gff=expand(os.path.join(workspace,
        "2.Annotation/ref/{ref}/{ref}.gff"), ref=REF),
        strain_gff=os.path.join(workspace,
        "2.Annotation/{group}_{strain}/{group}_{strain}.gff"),
    output:
        accessory_header=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory.header.embl"),
        accessory_tab=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory.tab"),
        accessory_binary_genes=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory_binary_genes.fa"),
        accessory_graph=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/accessory_graph.dot"),
        blast_identity_frequency=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/blast_identity_frequency.Rtab"),
        clustered_proteins=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/clustered_proteins"),
        core_accessory_header=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/core_accessory.header.embl"),
        core_accessory=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/core_accessory.tab"),
        core_accessory_graph=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/core_accessory_graph.dot"),
        core_alignment_header=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/core_alignment_header.embl"),
        core_gene_alignment=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/core_gene_alignment.aln"),
        gene_presence_absence_rtab=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/gene_presence_absence.Rtab"),
        gene_presence_absence_csv=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/gene_presence_absence.csv"),
        number_of_conserved_genes=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_conserved_genes.Rtab"),
        number_of_genes_in_pan_genome=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_genes_in_pan_genome.Rtab"),
        number_of_new_genes=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_new_genes.Rtab"),
        number_of_unique_genes=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/number_of_unique_genes.Rtab"),
        pan_genome_reference=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/pan_genome_reference.fa"),
        summary_statistics=os.path.join(workspace,
        "3.RoaryRef/{group}_{strain}_ref_pairwise/summary_statistics.txt"),
    message:
        "Run Roary for Strain-to-Ref 1:1 pairwise Pangenome analysis"
    params:
        threads=config["roary_pairwise_threads"],
        out_dir=os.path.join(workspace, "3.RoaryRef/{group}_{strain}_ref_pairwise/"),
        pident=config["roary_pairwise_pident"],
    conda:
        "env_roary"
    shell:
        """
        cd {params.out_dir}
        roary \
            -e \
            -p {params.threads} \
            -i {params.pident} \
            {input.ref_gff} {input.strain_gff}
        """


# Rule to move GFF files for within-group Roary 
rule move_gff_files:
    input:
        strain_gff=expand(os.path.join(workspace,
            "2.Annotation/{group}_{strain}/{group}_{strain}.gff"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
    output:
        tmp_dir=directory(os.path.join(workspace, "4.RoaryGroup/{group}/")),
    message:
        "Move GFF files for within-group Roary "
    params:
        workspace=config["workspace"],
        script=config["move_gff_script"],
        group_info=config["group_info"],
        tmp_dir=directory(os.path.join(workspace, "4.RoaryGroup/")),
    conda:
        "env_emapper"
    shell:
        """
        python workflow/scripts/shellmake.py \
            --group_yml config/{params.group_info} \
            --save_path {params.tmp_dir} \
            --script workflow/scripts/{params.script} \
            --workspace {params.workspace}

        bash workflow/scripts/{params.script}
        """


# Rule to run Roary for STRAINS within GROUPS
rule roary_within_group:
    input:
        # tmp_dir=os.path.join(workspace, "4.RoaryGroup/{group}/"),
        tmp_dir=rules.move_gff_files.output.tmp_dir,
    output:
        out_dir=directory(os.path.join(workspace, "4.RoaryGroup.tmp/{group}/")),
    params:
        threads=config["roary_within_group_threads"],
        pident=config["roary_within_group_pident"],
    message:
        "Run Roary for strains within groups"
    conda:
        "env_roary"
    shell:
        """
        mkdir {output.out_dir}
        cd {input.tmp_dir}

        roary \
            -e \
            -p {params.threads} \
            -i {params.pident} \
            *.gff
        
        mv *.gff {output.out_dir}
        """


#######################################################################################
# Rule to make six gene lists from Roary gene_presence_absence.csv
# core_all, core_hypo, core_nonhypo, shells_all, shells_hypo, shells_nonhypo
rule gene_list_maker:
    input:
        dummy=rules.roary_within_group.output.out_dir,
    output:
        out_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/gene_lists/")),
    message:
        "Make six gene lists from Roary gene_presence_absence.csv"
    params:
        tmp_dir=directory(os.path.join(workspace, "4.RoaryGroup/{group}/")),
    conda:
        "env_emapper"
    shell:
        """
        mkdir -p {output.out_dir}

        python workflow/scripts/gene_list_maker.py \
            --csv {params.tmp_dir}/gene_presence_absence.csv \
            --save_path {output.out_dir}/
        """


# Rule to copy FAA files from Prokka output to GROUP directories
rule move_faa_files:
    input:
        strain_faa=expand(os.path.join(workspace,
            "2.Annotation/{group}_{strain}/{group}_{strain}.faa"),
            zip,
            group=GROUP,
            strain=STRAIN,
        ),
    output:
        faa_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/faas/")),
    message:
        "Copy FAA files from Prokka output to group directories"
    params:
        workspace=config["workspace"],
        script=config["move_faa_script"],
        group_info=config["group_info"],
        group_dir=os.path.join(workspace, "5.AllGroups/"),
    conda:
        "env_emapper"
    shell:
        """
        python workflow/scripts/shellmake2.py \
            --group_yml config/{params.group_info} \
            --save_path {params.group_dir} \
            --script workflow/scripts/{params.script} \
            --workspace {params.workspace}

        bash workflow/scripts/{params.script}
        """


# Rule to curate six FASTA files based on gene lists created
rule fasta_curation:
    input:
        faa_dir=rules.move_faa_files.output.faa_dir,
        text_path=rules.gene_list_maker.output.out_dir,
    output:
        fasta_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/FASTA/")),
    message:
        "Curate six FASTA files based on gene lists created"
    params:
        tmp_dir=directory(os.path.join(workspace, "4.RoaryGroup/{group}/")),
        statistics=directory(os.path.join(workspace, "5.AllGroups/{group}/statistics.txt")),
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.fasta_dir}
        python workflow/scripts/core_shell_statistics_v2.py \
            --faa_path {input.faa_dir}/ \
            --text_path {input.text_path}/ \
            --save_path {output.fasta_dir}/ \
            --gpa {params.tmp_dir}/gene_presence_absence.csv \
            --summary {params.tmp_dir}/summary_statistics.txt > {params.statistics}
        """


# Rule to run COG analysis by running Egg-NOG-mapper for each group's core nonhypo FASTAs
rule cog_analysis_core:
    input:
        # Inputs are curated FASTAs from rule fasta_curation
        fasta_dir=rules.fasta_curation.output.fasta_dir,
    output:
        emapper_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/emapper-core/")),
    message:
        "Run COG analysis by running Egg-NOG-mapper for each group's core nonhypo FASTAs"
    params:
        # Parameters for Egg-NOG-mapper
        pident=config["emapper_core_pident"],
        evalue=config["emapper_core_evalue"],
        score=config["emapper_core_score"],
        query_cover=config["emapper_core_query_cover"],
        subject_cover=config["emapper_core_subject_cover"],
        output="{group}-core",
        cpu=config["emapper_core_cpu"],
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.emapper_dir}
        emapper.py \
            -i {input.fasta_dir}/core_nonhypo.fasta \
            --pident {params.pident} \
            --evalue {params.evalue} \
            --score {params.score} \
            --query_cover {params.query_cover} \
            --subject_cover {params.subject_cover} \
            --output {params.output} \
            --output_dir {output.emapper_dir} \
            --cpu {params.cpu}
        """


# Rule to run COG analysis by running Egg-NOG-mapper for each group's shell nonhypo FASTAs
rule cog_analysis_shells:
    input:
        # Inputs are curated FASTAs from rule fasta_curation
        fasta_dir=rules.fasta_curation.output.fasta_dir,
    output:
        emapper_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/emapper-shells/")),
    message:
        "Run COG analysis by running Egg-NOG-mapper for each group's shell nonhypo FASTAs"
    params:
        # Parameters for Egg-NOG-mapper
        pident=config["emapper_shells_pident"],
        evalue=config["emapper_shells_evalue"],
        score=config["emapper_shells_score"],
        query_cover=config["emapper_shells_query_cover"],
        subject_cover=config["emapper_shells_subject_cover"],
        output="{group}-shells",
        cpu=config["emapper_shells_cpu"],
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.emapper_dir}
        emapper.py \
            -i {input.fasta_dir}/shells_nonhypo.fasta \
            --pident {params.pident} \
            --evalue {params.evalue} \
            --score {params.score} \
            --query_cover {params.query_cover} \
            --subject_cover {params.subject_cover} \
            --output {params.output} \
            --output_dir {output.emapper_dir} \
            --cpu {params.cpu}
        """


# Rule to create nested pie chart for each group's core/shell COG analysis results
rule cog_visualization:
    input:
        # Inputs are annotated tsv files generated from Egg-NOG-mapper
        core_emapper_dir=rules.cog_analysis_core.output.emapper_dir,
        shells_emapper_dir=rules.cog_analysis_shells.output.emapper_dir,
    output:
        plot_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/cog_plots/")),
    message:
        "Create nested pie chart for each group's core/shell COG analysis results"
    params:
        core_hypo_path=os.path.join(workspace,
        "5.AllGroups/{group}/FASTA/core_hypo.fasta"),
        shells_hypo_path=os.path.join(workspace,
        "5.AllGroups/{group}/FASTA/shells_hypo.fasta"),
        group_name="{group}",
        types_core=config["types_core"],
        types_shells=config["types_shells"],
        core_statistics=os.path.join(workspace,
        "5.AllGroups/{group}/emapper-core/statistics.txt"),
        shells_statistics=os.path.join(workspace,
        "5.AllGroups/{group}/emapper-shells/statistics.txt"),
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.plot_dir}

        python workflow/scripts/cog_analysis_nested_plus.py \
            --tsv_file {input.core_emapper_dir}/*.annotations \
            --hypo_path {params.core_hypo_path} \
            --group_name {params.group_name} \
            --types {params.types_core} \
            --save_path {output.plot_dir}/ > {params.core_statistics}

        python workflow/scripts/cog_analysis_nested_plus.py \
            --tsv_file {input.shells_emapper_dir}/*.annotations \
            --hypo_path {params.shells_hypo_path} \
            --group_name {params.group_name} \
            --types {params.types_shells} \
            --save_path {output.plot_dir}/ > {params.shells_statistics}
        """


# Rule to create Roary plots - slightly modified version of roary_plots.py by Marco Galardini
# https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/README.md
rule roary_visualization:
    input:
        dummy=rules.roary_within_group.output.out_dir,
    output:
        plot_dir=directory(os.path.join(workspace, "5.AllGroups/{group}/roary_plots/")),
    message:
        "Create Roary plots"
    params:
        tmp_dir=directory(os.path.join(workspace, "4.RoaryGroup/{group}/")),
        group="{group}",
    conda:
        "env_emapper"
    shell:
        """
        mkdir {output.plot_dir}
        FastTree \
            -nt \
            -gtr {params.tmp_dir}/core_gene_alignment.aln > {params.group}_tree.newick

        python workflow/scripts/roary_plots.py \
            --labels {params.group}_tree.newick {params.tmp_dir}gene_presence_absence.csv \
            --save_path {output.plot_dir}

        mv {params.group}_tree.newick {output.plot_dir}
        """

        
# Rule to run Abricate to screen strain virulence factors
rule abricate_strain:
    input:
        # Input FNAs are annotated gene nucleotide file
        strain_fna=os.path.join(workspace, "2.Annotation/{group}_{strain}/{group}_{strain}.fna"),
    output:
        out_dir=directory(os.path.join(workspace, "6.Abricate/{group}_{strain}/")),
    message:
        "Run Abricate to screen strain virulence factors"
    params:
        minid=config["abricate_minid"],
        mincov=config["abricate_mincov"],
    conda:
        "env_abricate"
    shell:
        """
        mkdir {output.out_dir}

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db ncbi {input.strain_fna} > {output.out_dir}/ncbi.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db argannot {input.strain_fna} > {output.out_dir}/argannot.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db card {input.strain_fna} > {output.out_dir}/card.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db megares {input.strain_fna} > {output.out_dir}/megares.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db resfinder {input.strain_fna} > {output.out_dir}/resfinder.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db vfdb {input.strain_fna} > {output.out_dir}/vfdb.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db plasmidfinder {input.strain_fna} > {output.out_dir}/plasmidfinder.tab
        """


# Rule to run Abricate to screen reference virulence factors
rule abricate_ref:
    input:
        # Input FNAs are annotated gene nucleotide file
        ref_fna=os.path.join(workspace, "2.Annotation/ref/{ref}/{ref}.fna"),
    output:
        out_dir=directory(os.path.join(workspace, "6.Abricate/ref/{ref}/")),
    message:
        "Run Abricate to screen reference virulence factors"
    params:
        minid=config["abricate_minid"],
        mincov=config["abricate_mincov"],
    conda:
        "env_abricate"
    shell:
        """
        mkdir {output.out_dir}

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db ncbi {input.ref_fna} > {output.out_dir}/ncbi.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db argannot {input.ref_fna} > {output.out_dir}/argannot.tab
 
        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db card {input.ref_fna} > {output.out_dir}/card.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db megares {input.ref_fna} > {output.out_dir}/megares.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db resfinder {input.ref_fna} > {output.out_dir}/resfinder.tab

        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db vfdb {input.ref_fna} > {output.out_dir}/vfdb.tab
        
        abricate \
            --minid {params.minid} \
            --mincov {params.mincov} \
            --nopath \
            --db plasmidfinder {input.ref_fna} > {output.out_dir}/plasmidfinder.tab
        """