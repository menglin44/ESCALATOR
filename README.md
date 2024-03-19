# ESCALATOR
### polyg**E**nic **SC**ore h**A**rmonization and calcu**LAT**i**O**n sc**R**ipts
(Previous a.k.a. polyg**E**nic **S**core **CA**lcu**LAT**ion On eu**R**eka)

<img src="https://github.com/menglin44/ESCALATOR/assets/16557724/b5d68aa1-e18d-4a26-bd7d-751eace24011" width=20% height=20%>

This is a pipeline for harmonizing and calculating polygenic scores on a genetic dataset. 

* The original pipeline, specifically developed and to run on Eureka cloud platform for Colorado Center for Personalized Medicine, is archived under the subfolder [eureka_cloud_version](eureka_cloud_version)

* For general use on extensive platforms, we provide a containerized copy. Due to the file size of the container image, we do not host the image in this github, instead we offer

   - (1) the scripts and def file to build your container, with [instructions](https://github.com/MatthewFisher126/ESCALATOR?tab=readme-ov-file#using-the-container) from Matthew Fisher (matthew.j.fisher@cuanschutz.edu),

   - (2) or as an alternate, directly download the ready-to-use container image momentarily hosted [here](https://olucdenver-my.sharepoint.com/:u:/g/personal/meng_lin_cuanschutz_edu/ETTottyQgt5Akp3LkiORfFkBvmfnutRTTSHXQ3nlIAPIhg?e=VtdkJe).
 
(Please refer to the [separate folk](https://github.com/MatthewFisher126/ESCALATOR) for Matt's explanations and changes of generalized scripts and instructions of building the container)


### Overview

The pipeline takes care of build lifting, strand flipping, allele code mismatching etc. between external weights of PRS and the target genetic data, and calculation of final scores. 

<img src="https://github.com/menglin44/ESCALATOR/assets/16557724/4fa4ccf2-6f47-4a2c-a1bd-d96cf106b41c" width=70% height=70%>

### Usage

In light of using the containerized version, ESCALATOR can be run as 

```bash
singularity exec escalator-v1.sif masterPRS_format_v2_freeze3.sh [reformatting script designed (1, 2, or F)] \
[input directory (where weight file is)] \
[weight input filename] \
[output directory] \
[trait name (trait_PGSxxx)] \
[pfile directory] \
[pfile prefix name - ex: chr22_freeze3_dosages_PAIR.pgen = freeze3_dosages_PAIR]
```

Detailed explanations for logistics, along with usage examples, are described in the [vignette](escalator_container/ESCALATOR_container_readme.pdf) .







### Contact
Meng Lin (meng.lin@cuanschutz.edu) or Matthew Fisher (matthew.j.fisher@cuanschutz.edu)






