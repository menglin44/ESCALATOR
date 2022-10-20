# ESCALATOR
### polyg**E**nic **S**core **CA**lcu**LAT**ion On eu**R**eka 

This README file provides basic information re background and versions of the pipeline. Please see [prs_pipeline_readme.pdf](prs_pipeline_readme.pdf) for explanations and example usages.

The current version of directory hosts all scripts of the pipeline being used for quality control and calculation of polygenic scores of external weights for CCPM Biobank, which is imputed against TOPMed references and hosted/run on Google cloud based Eureka HPC. 

Briefly, the pipeline bridging takes care of build lifting, strand flipping, allele code mismatching etc. between external weights of PRS and the target genetic data, and calculation of final scores. 


This directory is temporary and being improved, and a more generic version of pipeline for customizable genetic data and non-cloud based environment will be added.

meng.lin@cuanschutz.edu
