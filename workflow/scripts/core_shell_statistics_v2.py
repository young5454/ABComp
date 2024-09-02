import os
import argparse
from Bio import SeqIO
import csv
from fasta_maker_v2 import fasta_maker_text, fasta_maker_text_v2

# 083124 : Revised version of core_shell_statistics.py

def core_all_fasta(faa_file, text_path, save_path):
    """
    Curates core_all.fasta from core_all.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    core_all_txt = os.path.join(text_path, 'core_all.txt')
    core_all_fasta = os.path.join(save_path, 'core_all.fasta')

    # Parse the faa_file and text_path to get query IDs
    with open(faa_file, 'r') as faa:
        faa_id_list = [record.id for record in SeqIO.parse(faa, "fasta")]

    with open(core_all_txt, 'r') as txt:
        core_all_id_list = txt.read().splitlines()

    # Save CDS core_all genes as FASTA
    cds_all_core_count = fasta_maker_text(fasta_original=faa_file,
                                          ids_text=core_all_txt,
                                          save_path=core_all_fasta)
    
    # These are core genes that are not protein-coding
    non_cds_core_count = 0
    no_entry_cores = []
    for query in core_all_id_list:
        if query not in faa_id_list:
            no_entry_cores.append(query)
            non_cds_core_count += 1
    
    return cds_all_core_count, non_cds_core_count, no_entry_cores


def core_nonhypo_fasta(faa_file, text_path, save_path):
    """
    Curates core_nonhypo.fasta from core_nonhypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    core_nonhypo_txt = os.path.join(text_path, 'core_nonhypo.txt')
    core_nonhypo_fasta = os.path.join(save_path, 'core_nonhypo.fasta')

    cds_core_non_hypo_count = fasta_maker_text(fasta_original=faa_file, 
                                               ids_text=core_nonhypo_txt, 
                                               save_path=core_nonhypo_fasta)

    return cds_core_non_hypo_count


def core_hypo_fasta(faa_file, text_path, save_path):
    """
    Curates core_hypo.fasta from core_hypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    core_hypo_txt = os.path.join(text_path, 'core_hypo.txt')
    core_hypo_fasta = os.path.join(save_path, 'core_hypo.fasta')

    cds_core_hypo_count = fasta_maker_text(fasta_original=faa_file, 
                                           ids_text=core_hypo_txt, 
                                           save_path=core_hypo_fasta)
    
    return cds_core_hypo_count


def shells_all_fasta(all_faa_files, text_path, save_path):
    """
    Curates shells_all.fasta from shells_all.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    shells_all_txt = os.path.join(text_path, 'shells_all.txt')
    shells_all_fasta = os.path.join(save_path, 'shells_all.fasta')

    # Parse the text_path to get query IDs
    with open(shells_all_txt, 'r') as txt:
        shells_all_id_list = txt.read().splitlines()

    # Loop through all FAA files to get all query IDs
    faa_id_list = []
    for faa_file in all_faa_files:
        with open(faa_file, 'r') as faa:
            for record in SeqIO.parse(faa, "fasta"):
                faa_id_list.append(record.id)
    
    # Loop through all FAA files and get matching sequences
    cds_all_shells_count = 0
    for faa_file in all_faa_files:
        # Save CDS shell genes as FASTA
        c = fasta_maker_text_v2(fasta_original=faa_file,
                                ids_text=shells_all_txt,
                                save_path=shells_all_fasta)
        cds_all_shells_count += c
        
    # These are shell genes that are not protein-coding
    non_cds_shells_count = 0
    no_entry_shells = []
    for query in shells_all_id_list:
        if query not in faa_id_list:
            no_entry_shells.append(query)
            non_cds_shells_count += 1
    
    return cds_all_shells_count, non_cds_shells_count, no_entry_shells
    

def shells_nonhypo_fasta(all_faa_files, text_path, save_path):
    """
    Curates shells_nonhypo.fasta from shells_nonhypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    shells_nonhypo_txt = os.path.join(text_path, 'shells_nonhypo.txt')
    shells_nonhypo_fasta = os.path.join(save_path, 'shells_nonhypo.fasta')

    # Loop through all FAA files and get matching sequences
    cds_shells_non_hypo_count = 0
    for faa_file in all_faa_files:
        # Save CDS shell genes as FASTA
        c = fasta_maker_text_v2(fasta_original=faa_file,
                                ids_text=shells_nonhypo_txt,
                                save_path=shells_nonhypo_fasta)
        cds_shells_non_hypo_count += c
        
    return cds_shells_non_hypo_count


def shells_hypo_fasta(all_faa_files, text_path, save_path):
    """
    Curates shells_hypo.fasta from shells_hypo.txt gene list
    Make sure to create proper gene list with gene_list_maker.py
    """
    shells_hypo_txt = os.path.join(text_path, 'shells_hypo.txt')
    shells_hypo_fasta = os.path.join(save_path, 'shells_hypo.fasta')
    
    # Loop through all FAA files and get matching sequences
    cds_shells_hypo_count = 0
    for faa_file in all_faa_files:
        # Save CDS shell genes as FASTA
        c = fasta_maker_text_v2(fasta_original=faa_file,
                                ids_text=shells_hypo_txt,
                                save_path=shells_hypo_fasta)
        cds_shells_hypo_count += c
        
    return cds_shells_hypo_count


def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Curate FASTAs of CDS core, shells genes. Also print result statistics')
    parser.add_argument('--faa_path', required=True, help='Path to FAA file')
    parser.add_argument('--text_path', required=True, help='Path to gene list')
    parser.add_argument('--save_path', required=True, help='Path to save gene FASTAs')
    parser.add_argument('--gpa', required=True, help='Path to gene_presence_absence.csv')
    parser.add_argument('--summary', required=True, help='Path to Roary summary_statistics.txt')

    args = parser.parse_args()

    faa_path = args.faa_path
    text_path = args.text_path
    save_path = args.save_path
    gpa = args.gpa
    summary = args.summary

    # Open the CSV file and output files
    with open(gpa, mode='r', newline='') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = list(csv_reader)[0]
        strains = header[14:]
    
    # List of all FAA files
    all_faa_files = [os.path.join(faa_path, f"{name}.faa") for name in strains]

    # First strain becomes the representative strain - Core genes will be identified with the representative strain IDs
    faa_file = all_faa_files[0]

    # Function calls - Cores
    cds_all_core_count, non_cds_core_count, no_entry_cores = core_all_fasta(faa_file=faa_file,
                                                                            text_path=text_path,
                                                                            save_path=save_path)
    cds_core_non_hypo_count = core_nonhypo_fasta(faa_file=faa_file, 
                                                 text_path=text_path, 
                                                 save_path=save_path)
    cds_core_hypo_count = core_hypo_fasta(faa_file=faa_file, 
                                          text_path=text_path,
                                          save_path=save_path)
    
    # Function calls - Shells
    cds_all_shells_count, non_cds_shells_count, no_entry_shells = shells_all_fasta(all_faa_files=all_faa_files,
                                                                                   text_path=text_path,
                                                                                   save_path=save_path)
    cds_shells_non_hypo_count = shells_nonhypo_fasta(all_faa_files=all_faa_files, 
                                                     text_path=text_path, 
                                                     save_path=save_path)
    cds_shells_hypo_count = shells_hypo_fasta(all_faa_files=all_faa_files, 
                                              text_path=text_path, 
                                              save_path=save_path)
    
    # Open summary_statistics.txt and parse Roary core & shell genes
    with open(summary, 'r') as file:
        # Split the data by lines
        lines = file.read().splitlines()

        # Extract Roary core & shell genes
        roary_core_genes = int(lines[0].split()[-1])
        roary_shell_genes = int(lines[2].split()[-1])

    # Print result statistics
    print('+--------------------------------------------------------------------+')
    print('+                          RESULT STATISTICS                         +')
    print('+--------------------------------------------------------------------+')
    print('+--------------------------------------------------------------------+')
    print('+        Successfully wrote all core queries to result FASTAs        +')
    print('+--------------------------------------------------------------------+')
    print('The total number of ROARY core genes are :', roary_core_genes)
    print('>> Please refer to summary_statistics.txt for this info \n')
    print('The total number of CDS all core genes are :', cds_all_core_count)
    print('The total number of CDS non-hypo core genes are :', cds_core_non_hypo_count)
    print('The total number of CDS hypo core genes are :', cds_core_hypo_count)
    print('The total number of non-CDS core genes are :', non_cds_core_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS non-hypo core + # of CDS hypo core = # of CDS all core?')
    sum = cds_core_non_hypo_count + cds_core_hypo_count
    print('>>', sum == cds_all_core_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS all core + # of non-CDS core = # of ROARY core?')
    sumy = cds_all_core_count + non_cds_core_count
    print('>>', sumy == roary_core_genes)
    print('+--------------------------------------------------------------------+')

    print('The ids of non-CDS core genes are as follows:')
    for ids in no_entry_cores: print(ids)
    print('These are one of the followings: tRNA, rRNA, tmRNA, ncRNA, ...')
    print('>> Please refer to .ffn file for each gene info')
        
    print('+--------------------------------------------------------------------+')
    print('+       Successfully wrote all shell queries to result FASTAs        +')
    print('+--------------------------------------------------------------------+')
    print('The total number of ROARY shell genes are :', roary_shell_genes)
    print('>> Please refer to summary_statistics.txt for this info \n')
    print('The total number of CDS all shell genes are :', cds_all_shells_count)
    print('The total number of CDS non-hypo shell genes are :', cds_shells_non_hypo_count)
    print('The total number of CDS hypo shell genes are :', cds_shells_hypo_count)
    print('The total number of non-CDS shell genes are :', non_cds_shells_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS non-hypo shells + # of CDS hypo shells = # of CDS all shells?')
    sumyy = cds_shells_non_hypo_count + cds_shells_hypo_count
    print('>>', sumyy == cds_all_shells_count)
    print('+--------------------------------------------------------------------+')

    print('# of CDS all shells + # of non-CDS shells = # of ROARY shells?')
    sumyyy = cds_all_shells_count + non_cds_shells_count
    print('>>', sumyyy == roary_shell_genes)
    print('+--------------------------------------------------------------------+')

    print('The ids of non-CDS shell genes are as follows:')
    for ids in no_entry_shells: print(ids)
    print('These are one of the followings: tRNA, rRNA, tmRNA, ncRNA, ...')
    print('>> Please refer to .ffn file for each gene info')
    print('+--------------------------------------------------------------------+')


if __name__ == "__main__":
    main()