#!/bin/sh

# system input variables
reform=${1} # using reformatting script designed for pgs catalog cleanup; 1, 2, or F
indir=${2}
infile=${3} # input weight file
dest=${4} # destination google bucket to store output
trait=${5} # trait name, if need to be left to decide as a trait_PGSxxx format from the input, input "unknown"
pfile_dir=${6} # Directory where the pfiles are
pfile=${7} # Name of plink fileset prefix AFTER chr# [ex: chr22_freeze3_dosages_PAIR.pgen = freeze3_dosages_PAIR]
# script_path=${7} # Full path to scripts
# bin_path=${8} # Full path to bin



echo "Working on input file ${infile}."

# predefined path and files
# pdir=''
# pfile='freeze3_dosages_PAIR'
script_path='/usr/bin'
bin_path='/usr/bin/'



#make sure dir doesn't have a last "/"
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

cd "${dest}"

# move appropriate data into tmp space from google bucket
# gsutil -q cp ${indir}/${infile} .
# gsutil -q cp ${script_path}/organize_pgs_format_v2.py .
# gsutil -q cp ${script_path}/makebed_v2.py .
# gsutil -q cp ${bin_path}/liftOver .
# gsutil -q cp ${bin_path}/*.chain.gz .
# gsutil -q cp ${script_path}/getname_format_v2.py .

# gsutil -q cp ${script_path}/atcg_bed.py .
# gsutil -q cp ${script_path}/bed2weightchr.py .
# gsutil -q cp ${bin_path}/plink2_mar .

# gsutil -q cp ${script_path}/integrate.R .

# chmod 700 ./*




##########################################
#@@@ Step I Preprocessing weight file @@@#
##########################################

# substract trait and PGS number to be part of the prefix
if [ "${trait}" = "unknown" ]
then
    echo "-----Step 0: Getting pgs name-----" 
    trait=$(/bin/python3 ${script_path}/getname_format_v2.py "${infile}") 
    echo -e "-----Done Step 0-----\n\n" 
else
    echo -e "-----Skipping Step 0: Getting pgs name-----\n\n"
fi

log="${trait}_prs.log"
touch "${log}"


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
        # gsutil -q cp ${log} ${dest}/
        # cp ${log} ${dest}/
        exit 1
    else
        echo "Reformatting of the input file failed." 2>&1 | tee -a "${log}"
        # gsutil -q cp ${log} ${dest}/
        # cp ${log} ${dest}/
        exit 1
    fi
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
        # gsutil -q cp ${log} ${dest}/
        # cp ${log} ${dest}/
        exit 1
    fi
fi

# quick summary
nsnps=$(wc -l "${indir}"/"${newinfile}" | awk '{print $1}')
echo "Total number of input weight file: ${nsnps}." 2>&1 | tee -a "${log}"

# convert to BED format: chr#, pos-1, pos, risk allele, other allele, weight
/bin/python3 ${script_path}/makebed_v2.py "${indir}"/"${newinfile}" "${dest}"/"${newinfile}".bed

# liftOver to hg38
if [[ "${lift}" != "F" ]]
then
    echo "-----Step 2: Lifting over from ${lift} to hg38-----" 2>&1 | tee -a "${log}"
    ${bin_path}/liftOver "${dest}"/"${newinfile}".bed ${bin_path}/${lift}ToHg38.over.chain.gz "${dest}"/"${newinfile}"_hg38.bed "${dest}"/"${newinfile}"_unmatched.txt
    bed="${newinfile}_hg38.bed"
    UNM=$(wc -l "${dest}"/"${newinfile}"_unmatched.txt | cut -d " " -f 1) # double lined
    UNM2=$((UNM/2))
    echo "Discarded unmatched variants from liftOver: ${UNM2}" 2>&1 | tee -a "${log}"
    echo -e "-----Done Step 2: Lifting over from ${lift} to hg38-----\n\n" 2>&1 | tee -a "${log}"
else
    echo "Assuming the file is already on hg38. No liftOver." 2>&1 | tee -a "${log}"
    echo -e "-----Skipping Step 2: Lifting over from ${lift} to hg38-----\n\n" 2>&1 | tee -a "${log}"
    bed="${newinfile}".bed
fi

# unifying filenames for the next step to fetch
mv "${dest}"/"${bed}" "${dest}"/"${trait}"_hg38_allchr.bed

# keeping only autosomal variants and count
for i in {21..22}
do
    awk '{if ($1=="chr""'"$i"'") print}' "${dest}"/"${trait}"_hg38_allchr.bed >> "${dest}"/"${trait}"_hg38.bed
done

nauto=$(wc -l "${dest}"/"${trait}"_hg38.bed | awk '{print $1}')
nall=$(wc -l "${dest}"/"${trait}"_hg38_allchr.bed | awk '{print $1}')
nnauto=$((nall - nauto))
echo "Number of non-autosomal variants being discarded: $nnauto" 2>&1 | tee -a "${log}"

#############################################
#@@@      Step II Running PRS per chr     @@@#
#############################################

echo "-----Step 3: Preparing to calculate PRS per chr-----" 2>&1 | tee -a "${log}"
echo "Will find overlapped loci in pvar files, remove at/cg, allele-mismatching, and flip strand accordingly. " 2>&1 | tee -a "${log}"

touch score.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tUpdatedSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest}"/"${trait}"_hg38_noAtCg_flipped.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest}"/"${trait}"_hg38_noAtCg_mismatch.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > "${dest}"/"${trait}"_hg38_noAtCg_missing_in_pvar.list
echo -e 'CHR(hg38)\tBP(hg38)\tOriginalSNPID\tUpdatedSNPID\tUpdatedRiskAllele\tUpdatedRefAllele\tWeight' > "${dest}"/"${trait}"_hg38_noAtCg_cleaned_forRecord.list

n_atcg_all=0 # keep the counts of at/cg loci
for CHR in {1..22}
do

    # subset to this chromosome
    awk 'BEGIN {OFS="\t"}{if ($1=="chr""'"${CHR}"'") print $0}' "${dest}"/"${trait}"_hg38.bed > "${dest}"/chr"${CHR}"_"${trait}"_hg38.bed
    # skip if no variants present
    nline=$(wc -l "${dest}"/chr"${CHR}"_"${trait}"_hg38.bed | awk '{print $1}')
    if [[ ${nline} == "0" ]]
    then
        continue
    fi 
    # cp the pgen file
    #gsutil -qm cp gs://hdchpcprodtis1-staging/mlin/freeze2_pgen/chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.* .
    # gsutil -qm cp ${pdir}/chr${CHR}_${pfile}.* .
    # remove ambiguous loci
    /bin/python3 ${script_path}/atcg_bed.py "${dest}"/chr"${CHR}"_"${trait}"_hg38.bed "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg.bed # also output another bed chr${CHR}_${trait}_hg38_at_cg.list
    bed="chr${CHR}_${trait}_hg38_noAtCg.bed"
    nout=$(wc -l "${dest}"/"${bed}" | awk '{print $1}')
    n_atcg=$((nline - nout))
    n_atcg_all=$((n_atcg_all + n_atcg))
    
    
    # match against pvar file, flip strand, remove unmatched / tri-allelic codes etc., and make into a weight file for PLINK
    #/bin/python3 bed2weightchr.py F ${bed} chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.pvar chr${CHR}_${trait}_hg38_noAtCg_clean.weights 
    /bin/python3 ${script_path}/bed2weightchr.py F "${dest}"/"${bed}" "${pfile_dir}"/chr"${CHR}"_${pfile}.pvar "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean.weights 
    ## also output chr${CHR}_${trait}_hg38_noAtCg_flipped.list,
    ## chr${CHR}_${trait}_hg38_noAtCg_mismatch.list,
    ## chr${CHR}_${trait}_hg38_noAtCg_missing_in_pvar.list
    ## chr${CHR}_${trait}_hg38_noAtCg_cleaned_forRecord.list

    # remove duplicated SNPs with the same risk allele but differerent weights (could be an error from the original weight or liftOver)
    cut -f 1,2 "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean.weights  | sort -k 1 | uniq -d | cut -f 1 > "${dest}"/chr"${CHR}"_duplicated.snps
    grep -vf "${dest}"/chr"${CHR}"_duplicated.snps "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean.weights  > "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_clean_noDup.weights
    weight="chr${CHR}_${trait}_hg38_noAtCg_clean_noDup.weights"
     
    # again, skip if no variants present
    nline=$(wc -l "${dest}"/"${weight}" | awk '{print $1}')
    if [[ "${nline}" == "0" ]]
    then
        # write snp info (flipped, discarded) to an integrated list
        cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_flipped.list >> "${dest}"/"${trait}"_hg38_noAtCg_flipped.list
        cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_mismatch.list >> "${dest}"/"${trait}"_hg38_noAtCg_mismatch.list
        cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_missing_in_pvar.list >> "${dest}"/"${trait}"_hg38_noAtCg_missing_in_pvar.list
        cat "${dest}"/chr"${CHR}"_duplicated.snps >> "${dest}"/"${trait}"_duplicated.list
        cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_cleaned_forRecord.list >> "${dest}"/"${trait}"_hg38_noAtCg_cleaned_forRecord.list
        
        #rm chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.*
        # rm chr${CHR}_${pfile}.* #uncomment
        continue
    fi 
    
    # Calculated PRS using PLINK2
    #./plink2_mar --pfile chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated \
    #             --score ${weight} list-variants \
    #             --out chr${CHR}_${trait}_prs
    ${bin_path}/plink2_mar --pfile "${pfile_dir}"/chr"${CHR}"_${pfile} \
                 --score "${dest}"/"${weight}" list-variants \
                 --out "${dest}"/chr"${CHR}"_${trait}_prs    
    
    
    # double checking if all weight SNPs were utilized
    nsnp_w=$(wc -l "${dest}"/"${weight}" | cut -d " " -f 1)
    nsnp_used=$(wc -l "${dest}"/chr"${CHR}"_${trait}_prs.sscore.vars | cut -d " " -f 1)
    if [[ "${nsnp_w}" != "${nsnp_used}" ]]
    then
        echo "There still are weight SNPs being discarded in PLINK --score after being cleaned on CHR${CHR}. Can check the log file for trouble shooting." 2>&1 | tee -a "${log}"
        #exit 1
    fi
    
    # remove the pgen file for step 3
    #rm chr${CHR}_freeze2_merged_overlapped_sites_INFOupdated.*
    # rm chr${CHR}_${pfile}.* #uncomment
    
    # write score name to a list
    echo "${dest}"/chr"${CHR}"_${trait}_prs.sscore >> "${dest}"/score.list
    # write snp info (flipped, discarded) to an integrated list
    cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_flipped.list >> "${dest}"/"${trait}"_hg38_noAtCg_flipped.list
    cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_mismatch.list >> "${dest}"/"${trait}"_hg38_noAtCg_mismatch.list
    cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_missing_in_pvar.list >> "${dest}"/"${trait}"_hg38_noAtCg_missing_in_pvar.list
    cat "${dest}"/chr"${CHR}"_duplicated.snps >> "${dest}"/"${trait}"_duplicated.list
    cat "${dest}"/chr"${CHR}"_"${trait}"_hg38_noAtCg_cleaned_forRecord.list >> "${dest}"/"${trait}"_hg38_noAtCg_cleaned_forRecord.list


done
echo -e "-----Done Step 3: Preparing to calculate PRS per chr-----\n\n" 2>&1 | tee -a "${log}"

###############################
#@@    Step III Integrating   @@#
###############################

echo -e "-----Step 4: Integrating scores from chromosomes.-----" 2>&1 | tee -a "${log}"
nscore=$(wc -l "${dest}"/score.list | awk '{print $1}')
if [[ "${nscore}" == "0" ]]
then
    echo "No variants are left for PRS score."
    # gsutil -q cp ${log} ${dest}/
    # cp ${log} ${dest}/
    exit 1
fi


Rscript ${script_path}/integrate.R "${dest}"/score.list "${dest}"/"${trait}"_prs.sscore

# organize the log report on # of snps
nflip=$(wc -l "${dest}"/${trait}_hg38_noAtCg_flipped.list | awk '{print $1}')
nmis=$(wc -l "${dest}"/${trait}_hg38_noAtCg_mismatch.list | awk '{print $1}')
nmissing=$(wc -l "${dest}"/${trait}_hg38_noAtCg_missing_in_pvar.list | awk '{print $1}')
ndup=$(wc -l "${dest}"/${trait}_duplicated.list | awk '{print $1}')
nused=$(wc -l "${dest}"/${trait}_hg38_noAtCg_cleaned_forRecord.list | awk '{print $1}')
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
# gsutil -q cp ${trait}_prs.sscore ${dest}/
# cp ${trait}_prs.sscore ${dest}/
# gsutil -q cp ${trait}_hg38_noAtCg*.list ${dest}/recordfiles/
cp "${dest}"/"${trait}"_hg38_noAtCg*.list "${dest}"/recordfiles/
# gsutil -q cp ${log} ${dest}/
# cp ${log} ${dest}/
# gsutil -q cp ${trait}_duplicated.list ${dest}/recordfiles/
cp "${dest}"/"${trait}"_duplicated.list "${dest}/"recordfiles/

