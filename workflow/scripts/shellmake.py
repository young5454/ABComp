import os
import argparse
import yaml
# 072424: Update for directory name change

parser = argparse.ArgumentParser(description="Shell script generator for moving GFF files")
parser.add_argument('--group_yml', required=True, help="Path to groups-strains information yaml file")
parser.add_argument('--save_path', required=True, help="Path to save GFF files")
parser.add_argument('--script', required=False, default="move_gffs.sh", help="Name of the shell script. Default is move_gff.sh")
parser.add_argument('--workspace', required=True, help="Name of the workspace directory")
args = parser.parse_args()

# Specify the path to your YAML file
groups_original = args.group_yml
save_path = args.save_path
script = args.script
workspace = args.workspace

# Read the YAML file and load it into a Python dictionary
with open(groups_original, 'r') as stream:
    try:
        groups_info = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Some logs for your help
print()
print("This is the group-strain info you provided:")
print(groups_info)
print()
print("Now making shell scripts for moving gff files into group folders...")

f = open(script, "w")
new_command = ""
main_path = os.path.join(workspace, "2.Annotation/")

# Make working directories
f.write(f"mkdir {save_path}\n")

for group in groups_info.keys():
    group_folder = os.path.join(save_path, f"{group}/")
    f.write(f"mkdir {group_folder}\n")
    for strain in groups_info[group]:
        gff_path = os.path.join(main_path, f"{group}_{strain}", f"{group}_{strain}.gff")
        new_command = f"cp {gff_path} {group_folder}\n"
        f.write(new_command)
f.close()

print()
print("Done")