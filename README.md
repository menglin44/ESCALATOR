# ESCALATOR README

Feel free to fork ESCALATOR itself from this github page or you can make it into a container. The tool/container will work in a SLURM environment or locally but it is recommended to use in a SLURM environment if your data set is large.  

## Using ESCALATOR by forking: 

**Note: You will need to fork this repo and edit the masterPRS_format_v2_freeze3.sh. Uncomment script_path and bin_path on lines 11 and 12 and comment out the same variables on lines 21 and 22.**

You will also need to unzip the prs_pipeline_bin.tar.gz file to access liftover and plink. 

Once the above step is complete, simply run the wrapper with the required arguments. 

An example is below:
```
# Arguments
Bash masterPRS_format_v2_freeze3.sh [reformatting script designed (1, 2, or F)] \
[input directory (where weight file is)] \
[weight input filename] \
[output directory] \
[trait name (trait_PGSxxx)] \
[pfile directory] \
[pfile prefix name - ex: chr22_freeze3_dosages_PAIR.pgen = freeze3_dosages_PAIR]
[path to scripts] \
[path to prs_pipeline_bin (plink and liftover)]

# Example
bash masterPRS_format_v2_freeze3.sh 2 \
./test_escalator/ \
non_hla_escalator_input.txt \
./out_escalator/ \
GRS2_non_hla \
./freeze3_w_dosages \
freeze3_dosages_PAIR \
./ESCALATOR/eureka_cloud_version/scripts/ \
./ESCALATOR/eureka_cloud_version/bin/prs_pipeline_bin/
``` 

## Using the container:

You can download the escalator.def file and create the container on https://cloud.sylabs.io/. You will also need singularity installed locally or on your HPC. 

Running the container is almost the same as above but you will specify the container name, the wrapper script, and then all of the arguments. Simply run the below:

```
singularity exec escalator-v1.sif masterPRS_format_v2_freeze3.sh [reformatting script designed (1, 2, or F)] \
[input directory (where weight file is)] \
[weight input filename] \
[output directory] \
[trait name (trait_PGSxxx)] \
[pfile directory] \
[pfile prefix name - ex: chr22_freeze3_dosages_PAIR.pgen = freeze3_dosages_PAIR]
```

You can also run the below to get similar information on running the container:

```
singularity run-help <container_name>.sif
```
