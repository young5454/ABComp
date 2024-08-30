import pandas as pd
import argparse
from fasta_maker import fasta_maker_list
    

def pairwise_analysis(csv_file, query_name, subject_name, query_faa, subject_faa, 
                      evalue_threshold, pident_threshold, qcovs_threshold, 
                      save_path, output_table_name, query_output_name, subject_output_name):
    print('_')
          
  
def main():    
    # Argument parser
    parser = argparse.ArgumentParser(
        description=
        """
        Code to parse type-to-type protein pairwise comparison analysis
        This code assumes Blastp analysis is conducted between two protein sets, with a qcovs option
        >> qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
        Input : Blastp output CSV file
        Output : Table of presence and absence of query proteins & FASTA file of pass query, subjects
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--csv_file', required=True, help='Path to pairwise Blastp output CSV file')
    parser.add_argument('-qn', '--query_name', required=True, help='Name of the query protein set')
    parser.add_argument('-sn', '--subject_name', required=True, help='Name of the subject protein set')
    parser.add_argument('-qf', '--query_faa', required=True, help='Path to representative FAA/FASTA file of query proteins')
    parser.add_argument('-sf', '--subject_faa', required=True, help='Path to representative FAA/FASTA file of subject proteins')
    parser.add_argument('-p', '--pident', required=False, default=80, help='Threshold for evalue. Default value is 80')
    parser.add_argument('-e', '--evalue', required=False, default=0.001, help='Threshold for percent identity. Default value is 0.001')
    parser.add_argument('-c', '--qcovs', required=False, default=70, help='Threshold for qcovs. Default value is 80')
    # parser.add_argument('--save_query', required=False, default='./', help='Path to save query pass FASTA. Defualt is current directory')
    # parser.add_argument('--save_subject', required=False, default='./', help='Path to save subject pass FASTA. Defualt is current directory')
    parser.add_argument('-s', '--save_path', required=False, default='./', help='Path to save final data')
    parser.add_argument('-o', '--output_table_name', required=False, default='presence_absence.csv', help='File name of presence and absence table')
    parser.add_argument('-qo', '--query_output_name', required=False, default='pass_query.fasta', help='File name of query pass FASTA. Defualt is pass_query.fasta')
    parser.add_argument('-so', '--subject_output_name', required=False, default='pass_subject.fasta', help='File name of subject pass FASTA. Defualt is pass_subject.fasta')
    args = parser.parse_args()
    
    # Define argument variables
    csv_file = args.csv_file

    query_name = args.query_name
    subject_name = args.subject_name
    query_faa = args.query_faa
    subject_faa = args.subject_faa

    evalue_threshold = float(args.evalue)
    pident_threshold = float(args.pident)
    qcovs_threshold = float(args.qcovs)
    save_path = args.save_path
    output_table_name = args.output_table
    query_output_name = args.query_output
    subject_output_name = args.subject_output

    # Run 
    table, pass_query, pass_subject = pairwise_analysis(csv_file=csv_file,
                                                        query_name=query_name,
                                                        subject_name=subject_name,
                                                        query_faa=query_faa,
                                                        subject_faa=subject_faa,
                                                        evalue_threshold=evalue_threshold,
                                                        pident_threshold=pident_threshold,
                                                        qcovs_threshold=qcovs_threshold,
                                                        save_path=save_path,
                                                        output_table_name=output_table_name,
                                                        query_output_name=query_output_name,
                                                        subject_output_name=subject_output_name)

    # Concatenate categories to the first line of CSV file
    category = "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore\n"

    # Read the existing contents of the CSV file
    with open(csv_file_path, 'r') as file:
        first_line = file.readline()
        existing_content = first_line + file.read()
        
        if first_line == category:
            updated_content = existing_content
        else:
            updated_content = category + existing_content

        # Write the updated content back to the CSV file
        with open(csv_file_path, 'w') as file:
            file.write(updated_content)

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file_path, sep='\t')
    query = df['qseqid']
    subject = df['sseqid']
    evalue = df['evalue']
    identity = df['pident']
    
    # Count query/subject entries that pass the identity threshold & eval threshold
    # Currently support only percent identity and e-value
    above_threshold = []
    pass_queries = []
    pass_subjects = []
    
    for i in range(len(identity)):
        curr_identity = identity[i]
        curr_evalue = evalue[i]
        
        if (curr_identity >= pident_threshold) and (curr_evalue < evalue_threshold):
            pass_query = query[i]
            pass_subject = subject[i]
            pass_queries.append(pass_query)
            pass_subjects.append(pass_subject)

    # Remove redundant entries
    pass_query_set = set(pass_queries)
    pass_subject_set = set(pass_subjects)

    # Print results
    print('+--------------------------------------------------------------------+')

    print('The number of unique query entries that pass the set threshold:', len(pass_query_set))
    print('The number of unique subject entries that pass the set threshold:', len(pass_subject_set))
    
    # Make FASTAs of passed entries (pass FASTAs)
    # Pass FASTA - query
    query_save_file_name = save_query + output_query
    query_count = fasta_maker_list(query_faa_file, pass_query_set, query_save_file_name)
    print('Successfully created and saved query pass FASTA file')
    
    # Pass FASTA - subject
    subject_save_file_name = save_subject + output_subject
    subject_count = fasta_maker_list(subject_faa_file, pass_subject_set, subject_save_file_name) 
    print('Successfully created and saved subject pass FASTA file')

    print('+--------------------------------------------------------------------+')

    
if __name__ == "__main__":
    main()
