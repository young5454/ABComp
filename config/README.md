# Configuration file
ABComp requires two configuration files for running the pipeline. These yaml files can be found in the `config/` directory. 

`config.yml` is a default configuration setting for the overall Snakemake run. Make sure you specify the correct parameters and directory names of your preference.

`groups_original.yml` is a configuration file for the complete group-strain information of your clinical isolates. Below is an example yaml file of a 2-group, 5-strain setting :

```YAML
NONMDR:
    - B0112
    - C0234
    - C3455
MDR:
    - B0232
    - D0991
```