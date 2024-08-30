# Default for Snakemake
echo Creating conda environment for default Snakemake : env_default ...
conda env create --name env_default --file workflow/envs/env_default.yml 

# Trimmomatic
echo Creating conda environment for Trimmomatic : env_trimmomatic ...
conda env create --name env_trimmomatic --file workflow/envs/env_trimmomatic.yml 

# FastQC and MultiQC
echo Creating conda environment for FastQC and MultiQC : env_qc ...
conda env create --name env_qc --file workflow/envs/env_qc.yml 

# Polypolish
echo Creating conda environment for Polypolish : env_polypolish ...
conda env create --name env_polypolish --file workflow/envs/env_polypolish.yml 

# Busco
echo Creating conda environment for Busco : env_busco ...
conda env create --name env_busco --file workflow/envs/env_busco.yml 

# Quast
echo Creating conda environment for Quast : env_quast ...
conda env create --name env_quast --file workflow/envs/env_quast.yml 

# Prokka
echo Creating conda environment for Prokka : env_prokka ...
conda env create --name env_prokka --file workflow/envs/env_prokka.yml 

# Roary
echo Creating conda environment for Roary : env_roary ...
conda env create --name env_roary --file workflow/envs/env_roary.yml 

# Egg-NOG-mapper & Downstream analyses
echo Creating conda environment for Egg-NOG-mapper and more downstream analyses : env_emapper ...
conda env create --name env_emapper --file workflow/envs/env_emapper.yml 

# Abricate
echo Creating conda environment for Abricate : env_abricate ...
conda env create --name env_abricate --file workflow/envs/env_abricate.yml 