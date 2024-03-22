#!/bin/sh

## system input variables
reform=${1} # using reformatting script designed for pgs catalog cleanup; 1, 2, 3, or F
indir=${2}
infile=${3} # input weight file
dest=${4} # destination google bucket to store output
trait=${5} # trait name, if need to be left to decide as a trait_PGSxxx format from the input, input "unknown"
pfile_dir=${6} # Directory where the pfiles are
pfile=${7} # Name of plink fileset midfix AFTER "chr#_" [ex: chr22_freeze3_dosages_PAIR.pgen = freeze3_dosages_PAIR]
rm_amb=${8} # T, or F to remove ambiguous variants. Default is T.
# script_path=${7} # Full path to scripts
# bin_path=${8} # Full path to bin

# make sure dir doesn't have a last "/"
if [ "${dest: -1}" = "/" ]
then
    temp="${dest%?}"
    dest="${temp}"
fi

if [ "${indir: -1}" = "/" ]
then
    temp="${indir%?}"
    indir="${temp}"
fi

if [ ! -d "${dest}" ]
then
    mkdir "${dest}"
fi
#cd "${dest}"


# make a temporary output folder with timestamp for intermediate files
starttime=`date +%s`
mkdir "${dest}"/"temp_${starttime}"
dest2="${dest}"/"temp_${starttime}"

# create a log file
log="${dest}"/"${trait}_prs.log"
touch "${log}"


## report params
echo "Input weight file ${indir}/${infile}." 2>&1 | tee -a "${log}"
if [ "${reform}" = "F" ]
then 
    echo "Weight file is indicated to be already reformatted. Will skip reformatting. " 2>&1 | tee -a "${log}"
else
    echo "Weight file is indicated to be in the format of catalog version ${reform}." 2>&1 | tee -a "${log}"
fi
echo "Will read autosomal genetic data from ${pfile_dir}/chr*_${pfile}.p*" 2>&1 | tee -a "${log}"
echo "Will write all output to ${dest}." 2>&1 | tee -a "${log}"
if [ "${rm_amb}" = "F" ]
then
    echo "User has opted for keeping variants with ambiguous allele code (e.g. A/T, C/G)." 2>&1 | tee -a "${log}"
else
    echo "Variants with ambiguous allele code (e.g. A/T, C/G) will be removed." 2>&1 | tee -a "${log}"
    rm_amb="T"
fi
echo '....................................................................................................' 2>&1 | tee -a "${log}"
echo -e '\n\n' 2>&1 | tee -a "${log}"
## end of report params



# predefined path for scripts and bin 
# pdir=''
# pfile='freeze3_dosages_PAIR'
script_path='/usr/bin'
bin_path='/usr/bin/'




##########################################
#@@@ Step I Preprocessing weight file @@@#
##########################################

# extract trait or PGS number to be part of the output prefix
if [ "${trait}" = "unknown" ]
then
    echo "-----Step 0: Getting pgs name-----" 
    trait=$(/bin/python3 ${script_path}/getname_format_v2.py "${infile}") 
    echo -e "-----Done Step 0-----\n\n" 
else
    echo -e "-----Skipping Step 0: Getting pgs name-----\n\n"
fi



# reformat input file if need to
if [ "${reform}" != "F" ]
then
    echo "-----Step 1: Reformatting the input weight file-----" 2>&1 | tee -a "${log}"
    /bin/python3 "${script_path}"/organize_pgs_format_v2.py "${reform}" "${indir}"/"${infile}"
    echo -e "-----Done Step 1-----\n\n" 2>&1 | tee -a "${log}"
    #newname=$(echo ${infile} | sed 's/\.txt/_reformated\.txt/g' | sed 's/\.gz//g')
    newname38=$(echo "${infile}" | sed 's/\.txt/_hg38_reformated\.txt/g' | sed 's/\.gz//g')
    newname19=$(echo "${infile}" | sed 's/\.txt/_hg19_reformated\.txt/g' | sed 's/\.gz//g')
    newname18=$(echo "${infile}" | sed 's/\.txt/_hg18_reformated\.txt/g' | sed 's/\.gz//g')
    newname17=$(echo "${infile}" | sed 's/\.txt/_hg17_reformated\.txt/g' | sed 's/\.gz//g')
    newnameNA=$(echo "${infile}" | sed 's/\.txt/_buildNA_reformated\.txt/g' | sed 's/\.gz//g')
    if [ -f "${indir}"/"${newname38}" ]
    then
        newinfile="${newname38}"
        lift="F"
    elif [ -f "${indir}"/"${newname19}" ]
    then
        newinfile="${newname19}"
        lift="hg19"
    elif [ -f "${indir}"/"${newname18}" ]
    then
        newinfile="${newname18}"
        lift="hg18"
    elif [ -f "${indir}"/"${newname17}" ]
    then
        newinfile="${newname17}"
        lift="hg17"
    elif [ -f "${indir}"/"${newnameNA}" ]
    then
        echo "Didn't detect the genome build of the input file." 2>&1 | tee -a "${log}"
        # gsutil -q cp ${log} ${dest2}/
        # cp ${log} ${dest2}/
        exit 1
    else
        echo "Reformatting of the input file failed." 2>&1 | tee -a "${log}"
        # gsutil -q cp ${log} ${dest2}/
        # cp ${log} ${dest2}/
        exit 1
    fi
    mv "${indir}"/"${newinfile}" "${dest2}"/"${newinfile}"
else
    echo -e "-----Skipping Step 1: Reformatting the input weight file-----\n\n" 2>&1 | tee -a "${log}"
    if [[ "${indir}"/"${infile}" == *"hg38"* ]]
    then
        newinfile="${infile}"
        lift="F"
    elif [[ "${indir}"/"${infile}" == *"hg19"* ]]
    then
        newinfile="${infile}"
        lift="hg19"
    elif [[ "${indir}"/"${infile}" == *"hg18"* ]]
    then
        newinfile="${infile}"
        lift="hg18"
    elif [[ "${indir}"/"${infile}" == *"hg17"* ]]
    then
        newinfile="${infile}"
        lift="hg17"
    else
        echo "Unknown genome build of the input file." 2>&1 | tee -a "${log}"
        # gsutil -q cp ${log} ${dest2}/
        # cp ${log} ${dest2}/
        exit 1
    fi
    cp "${indir}"/"${newinfile}" "${dest2}"/"${newinfile}"
fi

# quick summary
#nsnps=$(wc -l "${indir}"/"${newinfile}" | awk '{print $1}')
nsnps=$(wc -l "${dest2}"/"${newinfile}" | awk '{print $1}')
echo "Total number of input weight file: ${nsnps}." 2>&1 | tee -a "${log}"

# convert to BED format: chr#, pos-1, pos, risk allele, other allele, weight
/bin/python3 ${script_path}/makebed_v2.py "${dest2}"/"${newinfile}" "${dest2}"/"${newinfile}".bed

# liftOver to hg38
if [[ "${lift}" != "F" ]]
then
    echo "-----Step 2: Lifting over from ${lift} to hg38-----" 2>&1 | tee -a "${log}"
    ${bin_path}/liftOver "${dest2}"/"${newinfile}".bed ${bin_path}/${lift}ToHg38.over.chain.gz "${dest2}"/"${newinfile}"_hg38.bed "${dest2}"/"${newinfile}"_unmatched.txt
    bed="${newinfile}_hg38.bed"
    UNM=$(wc -l "${dest2}"/"${newinfile}"_unmatched.txt | cut -d " " -f 1) # double lined
    UNM2=$((UNM/2))
    echo "Discarded unmatched variants from liftOver: ${UNM2}" 2>&1 | tee -a "${log}"
    echo -e "-----Done Step 2: Lifting over from ${lift} to hg38-----\n\n" 2>&1 | tee -a "${log}"
else
    echo "Assuming the file is already on hg38. No liftOver." 2>&1 | tee -a "${log}"
    echo -e "-----Skipping Step 2: Lifting over from ${lift} to hg38-----\n\n" 2>&1 | tee -a "${log}"
    bed="${newinfile}".bed
fi

# unifying filenames for the next step to fetch
mv "${dest2}"/"${bed}" "${dest2}"/"${trait}"_hg38_allchr.bed

# keeping only autosomal variants and count
for i in {1..22}
do
    awk '{if ($1=="chr""'"$i"'") print}' "${dest2}"/"${trait}"_hg38_allchr.bed >> "${dest2}"/"${trait}"_hg38.bed
done

nauto=$(wc -l "${dest2}"/"${trait}"_hg38.bed | awk '{print $1}')
nall=$(wc -l "${dest2}"/"${trait}"_hg38_allchr.bed | awk '{print $1}')
nnauto=$((nall - nauto))
echo "Number of non-autosomal variants being discarded: $nnauto" 2>&1 | tee -a "${log}"

#############################################
#@@@      Step II Running PRS per chr     @@@#
#############################################

echo "-----Step 3: Preparing to calculate PRS per chr-----" 2>&1 | tee -a "${log}"
echo "Will find overlapped loci in pvar files, remove at/cg, allele-mismatching, and flip strand accordingly. " 2>&1 | tee -a "${log}"

touch ${dest2}/score.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tUpdatedSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest2}"/"${trait}"_flipped.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest2}"/"${trait}"_mismatch.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest2}"/"${trait}"_missing_in_pvar.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tUpdatedSNPID\tUpdatedRiskAllele\tUpdatedRefAllele\tWeight' > "${dest2}"/"${trait}"_cleaned_forRecord.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest2}"/"${trait}"_at_cg.list

n_atcg_all=0 # keep the counts of at/cg loci
for CHR in {1..22}
do

    # subset weight bed file to the current working chromosome
    awk 'BEGIN {OFS="\t"}{if ($1=="chr""'"${CHR}"'") print $0}' "${dest2}"/"${trait}"_hg38.bed > "${dest2}"/chr"${CHR}"_"${trait}"_hg38.bed
    # skip if no variants present
    nline=$(wc -l "${dest2}"/chr"${CHR}"_"${trait}"_hg38.bed | awk '{print $1}')
    if [[ ${nline} == "0" ]]
    then
        continue
    fi 
    # cp the pgen file
    #gsutil -qm cp gs://hdchpcprodtis1-staging/mlin/freeze2_pgen/chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.* .
    # gsutil -qm cp ${pdir}/chr${CHR}_${pfile}.* .
    
    # remove ambiguous loci (or not)
    if [[ "${rm_amb}" == "T" ]]
    then
        /bin/python3 ${script_path}/atcg_bed.py "${dest2}"/chr"${CHR}"_"${trait}"_hg38.bed "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg.bed # also output another bed chr${CHR}_${trait}_hg38_at_cg.list
        bed="chr${CHR}_${trait}_hg38_noAtCg.bed"
        nout=$(wc -l "${dest2}"/"${bed}" | awk '{print $1}')
        n_atcg=$((nline - nout))
        n_atcg_all=$((n_atcg_all + n_atcg))
    else
        mv "${dest2}"/chr"${CHR}"_"${trait}"_hg38.bed "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg.bed # for keeping downstream file names consistent - "noAtCg" in the name doesnt mean amb var removed
        touch "${dest2}"/chr${CHR}_${trait}_hg38_at_cg.list # an empty file
        bed="chr${CHR}_${trait}_hg38_noAtCg.bed"
        nout=$(wc -l "${dest2}"/"${bed}" | awk '{print $1}')
        n_atcg_all=0
    fi
    
    # match against pvar file, flip strand, remove unmatched / tri-allelic codes etc., and make into a weight file for PLINK
    #/bin/python3 bed2weightchr.py F ${bed} chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.pvar chr${CHR}_${trait}_hg38_noAtCg_clean.weights 
    /bin/python3 ${script_path}/bed2weightchr.py F "${dest2}"/"${bed}" "${pfile_dir}"/chr"${CHR}"_${pfile}.pvar "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean.weights 
    ## also output chr${CHR}_${trait}_hg38_noAtCg_flipped.list,
    ## chr${CHR}_${trait}_hg38_noAtCg_mismatch.list,
    ## chr${CHR}_${trait}_hg38_noAtCg_missing_in_pvar.list
    ## chr${CHR}_${trait}_hg38_noAtCg_cleaned_forRecord.list

    # remove duplicated SNPs with the same risk allele but differerent weights (could be an error from the original weight or liftOver)
    cut -f 1,2 "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean.weights  | sort -k 1 | uniq -d | cut -f 1 > "${dest2}"/chr"${CHR}"_duplicated.snps
    grep -vf "${dest2}"/chr"${CHR}"_duplicated.snps "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean.weights  > "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean_noDup.weights
    weight="chr${CHR}_${trait}_hg38_noAtCg_clean_noDup.weights"

     
    # again, skip if no variants present
    nline=$(wc -l "${dest2}"/"${weight}" | awk '{print $1}')
    if [[ "${nline}" == "0" ]]
    then
        # write snp info (flipped, discarded) to an integrated list
        cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_flipped.list >> "${dest2}"/"${trait}"_flipped.list
        cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_mismatch.list >> "${dest2}"/"${trait}"_mismatch.list
        cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_missing_in_pvar.list >> "${dest2}"/"${trait}"_missing_in_pvar.list
        cat "${dest2}"/chr"${CHR}"_duplicated.snps >> "${dest2}"/"${trait}"_duplicated.list
        cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_cleaned_forRecord.list >> "${dest2}"/"${trait}"_cleaned_forRecord.list
        cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_at_cg.list >> "${dest2}"/"${trait}"_at_cg.list
        #rm chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.*
        # rm chr${CHR}_${pfile}.* #uncomment
        continue
    fi 
    
    # Calculated PRS using PLINK2
    #./plink2_mar --pfile chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated \
    #             --score ${weight} list-variants \
    #             --out chr${CHR}_${trait}_prs
    ${bin_path}/plink2_mar --pfile "${pfile_dir}"/chr"${CHR}"_${pfile} \
                 --score "${dest2}"/"${weight}" list-variants no-mean-imputation \
                 --out "${dest2}"/chr"${CHR}"_${trait}_prs    


    # double checking if all weight SNPs were utilized
    nsnp_w=$(wc -l "${dest2}"/"${weight}" | cut -d " " -f 1)
    nsnp_used=$(wc -l "${dest2}"/chr"${CHR}"_${trait}_prs.sscore.vars | cut -d " " -f 1)
    if [[ "${nsnp_w}" != "${nsnp_used}" ]]
    then
        echo "There still are weight SNPs being discarded in PLINK --score after being cleaned on CHR${CHR}. Can check the log file for trouble shooting." 2>&1 | tee -a "${log}"
        #exit 1
    fi
    
    # remove the pgen file for step 3
    #rm chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.*
    # rm chr${CHR}_${pfile}.* #uncomment
    
    # write score name to a list
    echo "${dest2}"/chr"${CHR}"_${trait}_prs.sscore >> "${dest2}"/score.list
    # write snp info (flipped, discarded) to an integrated list
      cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_flipped.list >> "${dest2}"/"${trait}"_flipped.list
      cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_mismatch.list >> "${dest2}"/"${trait}"_mismatch.list
      cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_missing_in_pvar.list >> "${dest2}"/"${trait}"_missing_in_pvar.list
      cat "${dest2}"/chr"${CHR}"_duplicated.snps >> "${dest2}"/"${trait}"_duplicated.list
      cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_cleaned_forRecord.list >> "${dest2}"/"${trait}"_cleaned_forRecord.list
      cat "${dest2}"/chr"${CHR}"_"${trait}"_hg38_at_cg.list >> "${dest2}"/"${trait}"_at_cg.list
done
echo -e "-----Done Step 3: Preparing to calculate PRS per chr-----\n\n" 2>&1 | tee -a "${log}"

###############################
#@@    Step III Integrating   @@#
###############################

echo -e "-----Step 4: Integrating scores from chromosomes.-----" 2>&1 | tee -a "${log}"
nscore=$(wc -l "${dest2}"/score.list | awk '{print $1}')
if [[ "${nscore}" == "0" ]]
then
    echo "No variants are left for PRS score." 2>&1 | tee -a "${log}"
    # gsutil -q cp ${log} ${dest}/
    # cp ${log} ${dest}/
    exit 1
fi


Rscript ${script_path}/integrate.R "${dest2}"/score.list "${dest}"/"${trait}"_prs.sscore

# organize the log report on # of snps
nflip=$(wc -l "${dest2}"/${trait}_flipped.list | awk '{print $1}')
nmis=$(wc -l "${dest2}"/${trait}_mismatch.list | awk '{print $1}')
nmissing=$(wc -l "${dest2}"/${trait}_missing_in_pvar.list | awk '{print $1}')
ndup=$(wc -l "${dest2}"/${trait}_duplicated.list | awk '{print $1}')
nused=$(wc -l "${dest2}"/${trait}_cleaned_forRecord.list | awk '{print $1}')
echo "Number of ambigous A/T, C/G loci that are removed: ${n_atcg_all}" | tee -a "${log}"
echo "Number of weight SNPs that are flipped to the other strand: $((nflip-1))" | tee -a "${log}"
echo "Number of SNPs with mismatched allele codes against pgen that are removed: $((nmis-1))" | tee -a "${log}"
echo "Number of SNPs in weight file not found in pgen that are removed: $((nmissing-1))" | tee -a "${log}"
echo "Number of duplicated entries with same position and alleles (and are removed): ${ndup}" | tee -a "${log}"
echo "Number of final variants used for PRS: $((nused-1))" | tee -a "${log}"

echo -e "-----Done Step 4: Integrating scores from chromosomes.-----\n\n" 2>&1 | tee -a "${log}"
echo -e "Copying files to destinations..." 2>&1 | tee -a "${log}"

# Make recordfiles directory
if [[ ! -d "${dest}"/recordfiles/ ]]
then
    mkdir "${dest}"/recordfiles/
fi

# copy the files back to gs://
# gsutil -q cp ${trait}_prs.sscore ${dest2}/
# cp ${trait}_prs.sscore ${dest2}/
# gsutil -q cp ${trait}_hg38_noAtCg*.list ${dest2}/recordfiles/
#cp "${dest2}"/"${trait}"_hg38_noAtCg*.list "${dest2}"/recordfiles/
# gsutil -q cp ${log} ${dest2}/
# cp ${log} ${dest2}/
# gsutil -q cp ${trait}_duplicated.list ${dest2}/recordfiles/
#cp "${dest2}"/"${trait}"_duplicated.list "${dest2}/"recordfiles/

# mv record files to recrodfiles/ subfolder, instead of copying them 0315
mv "${dest2}"/"${trait}"_*.list "${dest}"/recordfiles/


# clean up 0315
rm -r "${dest2}"
