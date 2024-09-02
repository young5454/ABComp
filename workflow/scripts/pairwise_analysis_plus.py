import os
import sys
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
from fasta_maker_v2 import fasta_maker_list
    

def pairwise_analysis(raw_blastp_file, query_seq, subject_seq,
                      pident_threshold, evalue_threshold, qcovs_threshold,
                      save_path, best_match_table, not_founds_table, pass_query_fasta, pass_subject_fasta, not_founds_fasta):
        # Blastp Header
        header = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs"
        header = header.split('\t')

        # Read raw blastp output file
        raw_blastp_file = pd.read_csv(raw_blastp_file, sep="\t", names=header)
        ## Print log
        print('+--------------------------------------------------------------------+')
        print('Total number of raw entries :', len(raw_blastp_file))

        # Filtering with thresholds
        # 1. Percent identity threshold
        fil_blastp_file = raw_blastp_file[raw_blastp_file['pident'] >= pident_threshold]
        ## Print log
        print('+--------------------------------------------------------------------+')
        print('Total number of entries with pident threshold applied :', len(fil_blastp_file))

        # 2. e-value threshold
        fil_blastp_file = fil_blastp_file[fil_blastp_file['evalue'] < evalue_threshold]
        print('+--------------------------------------------------------------------+')
        print('Total number of entries with e-value threshold applied :', len(fil_blastp_file))

        # 3. Query coverage threshold
        fil_blastp_file = fil_blastp_file[fil_blastp_file['qcovs'] >= qcovs_threshold]
        print('+--------------------------------------------------------------------+')
        print('Total number of entries with query coverage threshold applied :', len(fil_blastp_file))

        # Get Not-Found queries
        all_unique_queries = np.unique(raw_blastp_file['qseqid'])
        founds_list = np.unique(fil_blastp_file['qseqid'])
        not_founds_list = [prot for prot in all_unique_queries if prot not in founds_list]

        # Terminate if 0 founds
        if len(founds_list) == 0:
             print('Zero number of founds : Terminate...')
             sys.exit()

        # Prepare outputs & save them to save_path
        # 1. Best Matches table
        best_matches = fil_blastp_file.drop_duplicates(subset='qseqid', keep='first')
        best_matches.to_csv(os.path.join(save_path, best_match_table), sep="\t", index=False)

        # 2. Not-Founds table
        not_founds = raw_blastp_file[raw_blastp_file['qseqid'].isin(not_founds_list)]
        not_founds.to_csv(os.path.join(save_path, not_founds_table), sep="\t", index=False)

        # 3. Pass query FASTA
        pass_query_list = list(best_matches['qseqid'])
        fasta_maker_list(fasta_original=query_seq,
                         ids_list=pass_query_list,
                         save_path=os.path.join(save_path, pass_query_fasta))
        
        # 4. Pass subject FASTA
        pass_subject_list = list(best_matches['sseqid'])
        fasta_maker_list(fasta_original=subject_seq,
                         ids_list=pass_subject_list,
                         save_path=os.path.join(save_path, pass_subject_fasta))
        
        # 5. Not-Founds FASTA
        fasta_maker_list(fasta_original=query_seq,
                         ids_list=not_founds_list,
                         save_path=os.path.join(save_path, not_founds_fasta))

        print("All outputs saved successfully")
        print('+--------------------------------------------------------------------+')
          
  
def main():    
    # Argument parser
    parser = argparse.ArgumentParser(
        description=
        """
        ++ This code will ease comparing two protein sets using Blastp ++
        After a manual Blastp analysis is done, pairwise_analysis_plus.py will do the following :
        1. Parse raw Blastp output table with pident, evalue and qcovs threshold for determining presence of query entry within subject set
        2. The single best query-to-subject match (Best Matches) will be saved as a Blastp output table
        3. The best match queries and subjects will be saved as a FASTA file 
        4. Not-Found queries will be saved as a Blastp output table with all results
        5. Not-Found queries will be saved as a FASTA file

        This code assumes Blastp analysis is conducted between two protein sets, with a qcovs option
        >> qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # Inputs
    parser.add_argument('-b', '--raw_blastp_file', required=True, help='Path to raw Blastp output file of pairwise comparison between two protein sets')
    parser.add_argument('-qf', '--query_seq', required=True, help='Path to representative FAA/FASTA file of query proteins')
    parser.add_argument('-sf', '--subject_seq', required=True, help='Path to representative FAA/FASTA file of subject proteins')

    # Thresholds
    parser.add_argument('-p', '--pident', required=False, default=35, help='Threshold for percent identity. Default value is 35')
    parser.add_argument('-e', '--evalue', required=False, default=0.001, help='Threshold for evalue. Default value is 0.001')
    parser.add_argument('-c', '--qcovs', required=False, default=70, help='Threshold for qcovs. Default value is 70')

    # Outputs
    parser.add_argument('-s', '--save_path', required=False, default='./', help='Path to save all final data')
    parser.add_argument('-bmt', '--best_match_table', required=False, default='best_match.csv', help='File name of best match Blastp output table')
    parser.add_argument('-nft', '--not_founds_table', required=False, default='not_founds.csv', help='File name of not-founds table')
    parser.add_argument('-pqf', '--pass_query_fasta', required=False, default='pass_query.fasta', help='File name of best match query FASTA')
    parser.add_argument('-psf', '--pass_subject_fasta', required=False, default='pass_subject.fasta', help='File name of best match subject FASTA')
    parser.add_argument('-nff', '--not_founds_fasta', required=False, default='not_founds.fasta', help='File name of not-founds query FASTA')

    args = parser.parse_args()
    
    # Define argument variables
    raw_blastp_file = args.raw_blastp_file
    query_seq = args.query_seq
    subject_seq = args.subject_seq

    pident_threshold = float(args.pident)
    evalue_threshold = float(args.evalue)
    qcovs_threshold = float(args.qcovs)

    save_path = args.save_path
    best_match_table = args.best_match_table
    not_founds_table = args.not_founds_table
    pass_query_fasta = args.pass_query_fasta
    pass_subject_fasta = args.pass_subject_fasta
    not_founds_fasta = args.not_founds_fasta

    # Run 
    pairwise_analysis(raw_blastp_file=raw_blastp_file,
                      query_seq=query_seq,
                      subject_seq=subject_seq,
                      pident_threshold=pident_threshold,
                      evalue_threshold=evalue_threshold,
                      qcovs_threshold=qcovs_threshold,
                      save_path=save_path,
                      best_match_table=best_match_table,
                      not_founds_table=not_founds_table,
                      pass_query_fasta=pass_query_fasta,
                      pass_subject_fasta=pass_subject_fasta,
                      not_founds_fasta=not_founds_fasta)

    
if __name__ == "__main__":
    main()
