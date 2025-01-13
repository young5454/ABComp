import os
import yaml
from Bio import SeqIO
import argparse


def init_assembly_dir(groups_original_yml, base_dir):
    """
    Initialize the 0.Assembly directory structure using groups_original.yml file
    """
    with open(groups_original_yml, 'r') as file:
        data = yaml.safe_load(file)
    
    # Group-Strain directories
    for group, strains in data.items():
        for strain in strains:
            strain_dir = os.path.join(base_dir, f'{group}_{strain}')
            genome_dir = os.path.join(strain_dir, 'genome')
            reads_dir = os.path.join(strain_dir, 'reads')
            
            # Create directories
            os.makedirs(genome_dir, exist_ok=True)
            os.makedirs(reads_dir, exist_ok=True)

    # Reference directory
    ref_dir = os.path.join(base_dir, 'ref')
    ref_genome_dir = os.path.join(ref_dir, 'genome')
    os.makedirs(ref_genome_dir, exist_ok=True)
    

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Initialize the 0.Assembly directory structure using groups_original.yml file')
    parser.add_argument('-g', '--groups', required=True, help='Path to groups_original.yml file')
    parser.add_argument('-b', '--base', required=True, help='Path to 0.Assembly inside workspace directory')

    args = parser.parse_args()

    # Define argument variables
    groups = args.groups
    base = args.base

    # Run
    init_assembly_dir(groups_original_yml=groups,
                      base_dir=base)


if __name__ == "__main__":
    main()
