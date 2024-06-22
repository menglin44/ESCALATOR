# ESCALATOR

$${\color{olive}polyg**E**nic \ **SC**ore \ h**A**rmonization \ and \ calcu**LAT**i**O**n \ sc**R**ipts}$$

$${\color{grey}(Previously \ a.k.a. \ polyg**E**nic \ **S**core \ **CA**lcu**LAT**ion \ On \ eu**R**eka)}$$

<img src="https://github.com/menglin44/ESCALATOR/assets/16557724/b5d68aa1-e18d-4a26-bd7d-751eace24011" width=20% height=20%>

## Overview

This is a pipeline for harmonizing and calculating polygenic scores in a genetic dataset. It takes care of build lifting, strand flipping, allele code mismatching etc. between external weights of PRS and the target genetic data, and calculation of final scores. 

<img src="https://github.com/menglin44/ESCALATOR/assets/16557724/13a5464e-f111-425d-bdcb-caed6ab68ca5" width=70% height=70%>

## Download

There are two ways to use ESCALATOR -

* For general use on extensive platforms, we provide a singularity container of ESCALATOR. It requires that [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) is already installed on your platform. 

   - Download the singularity container image [escalator-v2.sif](https://olucdenver-my.sharepoint.com/:u:/g/personal/meng_lin_cuanschutz_edu/EQ8IM0p0itZHgKGqKge6JY0BVXAovZ66TpeV6waKr100DQ)
   
   - Alternatively, you can also easily build the container yourself using the original scripts hosted at this [repo](https://github.com/MatthewFisher126/ESCALATOR), and follow the instructions with the def file by Matthew Fisher (matthew.j.fisher@cuanschutz.edu).

   - If you prefer a docker version of container than singularity, contact matthew.j.fisher@cuanschutz.edu

* If you wish to use the original scripts, which uses a bash wrapper script to call other helpers: You can download the scripts from the script holder repo above, and will need to lightly modify the script by following the instructions [here](https://github.com/MatthewFisher126/ESCALATOR?tab=readme-ov-file#if-to-use-the-original-scripts-by-download-the-repo)

> [!WARNING]
> The deprecated original pipeline scripts, initially developed for Eureka cloud platform for Colorado Center for Personalized Medicine, are archived under the subfolder [eureka_hpc_deprecated](eureka_hpc_deprecated). These deprecated scripts are archived as a test prototype and **should not** be used.


## Usage
> [!Important]
> Complete explanations of how to use ESCALATOR can be found in the [vignette](escalator_container/ESCALATOR_container_readme.pdf)

In brief of the container's usage, ESCALATOR can be run as 

```
singularity exec escalator-v2.sif masterPRS_v4.sh [reformatting script designed (1, 2, or 3)] \
[input directory (where weight file is)] \
[weight input filename] \
[output directory] \
[trait name (trait_PGSxxx)] \
[pfile directory] \
[pfile prefix name - ex: chr22_freeze3_dosages_PAIR.pgen = freeze3_dosages_PAIR] \
[T or F - whether to remove ambiguous variants] \
[NA or filename - frequency file for PLINK to impute missing genotypes, can be NA to skip if sample size >50]
```

### Run with an example

> [!Tip]
> You can download and decompress the example data in [dosages](https://github.com/menglin44/ESCALATOR/blob/main/test%20data/example_dosages.tar.gz) or in [hardcall genotypes with missingness]([test_data/example_hardcalls.tar.gz](https://github.com/menglin44/ESCALATOR/blob/main/test%20data/example_hardcalls.tar.gz)), and pair with the any of the [weight file input](https://github.com/menglin44/ESCALATOR/tree/main/test%20data/example%20weight%20input) for a test run.

An example for running with all data available in a folder presumably called /home/test/, becomes

```bash
singularity exec escalator-v2.sif masterPRS_v4.sh 2 \
/home/test/ \
dummy_lv5.wgt.txt \
/home/test/ \
testrun \
/home/test/example_dosages/ \
dummy_dosage \
T \
NA
```
Here the command tells the ESCALATOR to remove any amibigous variants in the weight file.




### Contact
Meng Lin (meng.lin@cuanschutz.edu) or Matthew Fisher (matthew.j.fisher@cuanschutz.edu)






