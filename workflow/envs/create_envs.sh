# Default for Snakemake
echo Creating conda environment for default Snakemake: env_default ...
conda env create --name env_default --file env_default.yml 

# Polypolish
echo Creating conda environment for Polypolish: env_polypolish ...
conda env create --name env_polypolish --file env_polypolish.yml 

# Busco
echo Creating conda environment for Busco: env_busco ...
conda env create --name env_busco --file env_busco.yml 

# Quast
echo Creating conda environment for Quast: env_quast ...
conda env create --name env_quast --file env_quast.yml 

# Prokka
echo Creating conda environment for Prokka: env_prokka ...
conda env create --name env_prokka --file env_prokka.yml 

# Roary
echo Creating conda environment for Roary: env_roary ...
conda env create --name env_roary --file env_roary.yml 

# Egg-NOG-mapper & Downstream analyses
echo Creating conda environment for Egg-NOG-mapper and more downstream analyses: env_emapper ...
conda env create --name env_emapper --file env_emapper.yml 