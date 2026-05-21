
# Setup environment
export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2026/bin/:$PATH

# Go to work directory
cd $VSC_SCRATCH/snakemake_snpcalling 

# Link input data
ln -s /lustre1/project/stg_00079/teaching/I0U19a_data_2026/TLE66_N.fastq .
ln -s /lustre1/project/stg_00079/teaching/I0U19a_data_2026/TLE66_T.fastq .
ln -s /staging/leuven/stg_00079/teaching/hg38_9/chr9.fa .

# dry run
snakemake --snakefile /data/leuven/381/vsc38194/r1079125_Miaojun_Zhang/030_snakemake_snpcalling/Snakefile -n
# Run snakemake
using the `--snakefile` parameter pointing to original snakemake file

snakemake --snakefile /data/leuven/381/vsc38194/r1079125_Miaojun_Zhang/030_snakemake_snpcalling/Snakefile -c 4

snakemake -c1 --snakefile Snakefile

# search for the path to the original snakefile 
find $VSC_HOME -name "Snakefile" 2>/dev/null
find $VSC_DATA -name "Snakefile" 2>/dev/null

# connect to Mango
iron auth <user name> <zone name> <host>
iron auth r1029125 	gbiomed gbiomed.irods.icts.kuleuven.be

# relative paths
current directory (.), parent directory (..), your
home director (˜), someone else’s home directory (˜ma51k),
relative to current location (scs/staﬀ if you are at /home)

# git push 
check status: git status
git add file / git add .
git commit -m "message"
git push
