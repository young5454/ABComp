from Bio import SeqIO

# 083124 : Revised version of fasta_maker.py

def fasta_maker_list(fasta_original, ids_list, save_path):
    """
    Makes a new FASTA file by finding ids of ids_list in fasta_original
    Saves the new FASTA to save_path
    """
    result_sequences = []

    # Parse the input FASTA file using SeqIO
    with open(fasta_original, 'r') as faa:
        for record in SeqIO.parse(faa, "fasta"):
            # Check if the record's ID is in the provided list of IDs
            if record.id in ids_list:
                result_sequences.append(record)

    # Write the matched sequences to the output FASTA file
    with open(save_path, 'w') as result_file:
        SeqIO.write(result_sequences, result_file, "fasta")


def fasta_maker_text(fasta_original, ids_text, save_path):
    """
    Makes a new FASTA file by finding ids of ids_text in fasta_original
    The ID text file assumes one id per line (with newline)
    Saves the new FASTA to save_path
    """
    result_sequences = []

    # Parse the ids_text to get ids
    with open(ids_text, 'r') as txt:
        ids_list = txt.read().splitlines()

    count = 0
    # Parse the input FASTA file using SeqIO
    with open(fasta_original, 'r') as faa:
        for record in SeqIO.parse(faa, "fasta"):
            # Check if the record's ID is in the provided list of IDs
            if record.id in ids_list:
                result_sequences.append(record)
                count += 1

    # Write the matched sequences to the output FASTA file
    with open(save_path, 'w') as result_file:
        SeqIO.write(result_sequences, result_file, "fasta")
    
    return count


def fasta_maker_text_v2(fasta_original, ids_text, save_path):
    """
    Makes a new FASTA file by finding ids of ids_text in fasta_original
    The ID text file assumes one id per line (with newline)
    Append instead of overwriting
    """
    result_sequences = []

    # Parse the ids_text to get ids
    with open(ids_text, 'r') as txt:
        ids_list = txt.read().splitlines()

    count = 0
    # Parse the input FASTA file using SeqIO
    with open(fasta_original, 'r') as faa:
        for record in SeqIO.parse(faa, "fasta"):
            # Check if the record's ID is in the provided list of IDs
            if record.id in ids_list:
                result_sequences.append(record)
                count += 1

    # Write the matched sequences to the output FASTA file - Append instead of overwriting
    with open(save_path, 'a') as result_file:
        SeqIO.write(result_sequences, result_file, "fasta")
    
    return count