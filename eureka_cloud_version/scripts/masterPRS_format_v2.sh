#!/bin/sh

# default path to files (set on augrabies)
SPATH="/share/hennlab/lab_scripts/prs_pipeline_lemon"
REFPATH="/share/hennlab/reference/liftover_references"
BINPATH="/share/hennlab/progs"
BPATH="/share/hennlab/data/snp-array/UKB_arrays/DanaGil"
BFILE="SoAs_imp90"
BUILD="Hg19"
dest="."


# system input variables
reform=${1} # using reformatting script designed for pgs catalog cleanup; 1, 2, or F
indir=${2}
infile=${3} # input weight file
#dest=${4} # destination google bucket to store output
trait=${4} # trait name, if need to be left to decide as a trait_PGSxxx format from the input, input "unknown"


echo "Working on input file ${infile}."

#make sure dir doesn't have a last "/"
if [ ${dest: -1} = "/" ]
then
    temp=${dest%?}
    dest=${temp}
fi

if [ ${indir: -1} = "/" ]
then
    temp=${indir%?}
    indir=${temp}
fi

# move appropriate data into tmp space from google bucket
cp ${indir}/${infile} .
#cp ${SPATH}/organize_pgs_format_v2.py .
#cp ${SPATH}/makebed_v2.py .
#cp ${BINPATH}/liftOver .
#cp ${REFPATH}/*.chain.gz .
#cp ${SPATH}/getname_format_v2.py .

#cp ${SPATH}/atcg_bed.py .
#cp ${SPATH}/bed2weightchr.py .
#cp ${BINPATH}/plink2_mar .

#cp ${SPATH}/integrate.R .

chmod 755 ./*




##########################################
#@@@ Step I Preprocessing weight file @@@#
##########################################

# substract trait and PGS number to be part of the prefix
if [ ${trait} = "unknown" ]
then
    echo "-----Step 0: Getting pgs name-----" 
    trait=$(python2 ${SPATH}/getname_format_v2.py ${infile}) 
    echo -e "-----Done Step 0-----\n\n" 
else
    echo -e "-----Skipping Step 0: Getting pgs name-----\n\n"
fi

log="${trait}_prs.log"
touch ${log}


# reformat input file if need to
if [ ${reform} != "F" ]
then
    echo "-----Step 1: Reformatting the input weight file-----" 2>&1 | tee -a ${log}
    python2 ${SPATH}/organize_pgs_format_v2.py ${reform} ${infile}
    echo -e "-----Done Step 1-----\n\n" 2>&1 | tee -a ${log}
    #newname=$(echo ${infile} | sed 's/\.txt/_reformated\.txt/g' | sed 's/\.gz//g')
    newname38=$(echo ${infile} | sed 's/\.txt/_hg38_reformated\.txt/g' | sed 's/\.gz//g')
    newname19=$(echo ${infile} | sed 's/\.txt/_hg19_reformated\.txt/g' | sed 's/\.gz//g')
    newname18=$(echo ${infile} | sed 's/\.txt/_hg18_reformated\.txt/g' | sed 's/\.gz//g')
    newname17=$(echo ${infile} | sed 's/\.txt/_hg17_reformated\.txt/g' | sed 's/\.gz//g')
    newnameNA=$(echo ${infile} | sed 's/\.txt/_buildNA_reformated\.txt/g' | sed 's/\.gz//g')
    if [ -f ${newname38} ]
    then
        newinfile="${newname38}"
        lift="hg38"
    elif [ -f ${newname19} ]
    then
        newinfile="${newname19}"
        lift="hg19"
    elif [ -f ${newname18} ]
    then
        newinfile="${newname18}"
        lift="hg18"
    elif [ -f ${newname17} ]
    then
        newinfile="${newname17}"
        lift="hg17"
    elif [ -f ${newnameNA} ]
    then
        echo "Didn't detect the genome build of the input file." 2>&1 | tee -a ${log}
        cp ${log} ${dest}/
        exit 1
    else
        echo "Reformatting of the input file failed." 2>&1 | tee -a ${log}
        cp ${log} ${dest}/
        exit 1
    fi
else
    echo -e "-----Skipping Step 1: Reformatting the input weight file-----\n\n" 2>&1 | tee -a ${log}
    if [[ "${infile}" == *"hg38"* ]]
    then
        newinfile=${infile}
        lift="F"
    elif [[ "${infile}" == *"hg19"* ]]
    then
        newinfile=${infile}
        lift="hg19"
    elif [[ "${infile}" == *"hg18"* ]]
    then
        newinfile=${infile}
        lift="hg18"
    elif [[ "${infile}" == *"hg17"* ]]
    then
        newinfile=${infile}
        lift="hg17"
    else
        echo "Unknown genome build of the input file." 2>&1 | tee -a ${log}
        cp ${log} ${dest}/
        exit 1
    fi
fi

# quick summary
nsnps=$(wc -l ${newinfile} | awk '{print $1}')
echo "Total number of input weight file: ${nsnps}." 2>&1 | tee -a ${log}

# convert to BED format: chr#, pos-1, pos, risk allele, other allele, weight
python2 ${SPATH}/makebed_v2.py ${newinfile} ${newinfile}.bed

# liftOver to hg38
temp_build=$(echo "$BUILD" | awk '{print tolower($0)}')
if [[ ${lift} != ${temp_build} ]]
then
    echo "-----Step 2: Lifting over from ${lift} to desired build-----" 2>&1 | tee -a ${log}
    ${BINPATH}/liftOver ${newinfile}.bed ${REFPATH}/${lift}To${BUILD}.over.chain.gz ${newinfile}_${BUILD}.bed ${newinfile}_unmatched.txt
    bed="${newinfile}_${BUILD}.bed"
    UNM=$(wc -l ${newinfile}_unmatched.txt | cut -d " " -f 1) # double lined
    UNM2=$((UNM/2))
    echo "Discarded unmatched variants from liftOver: ${UNM2}" 2>&1 | tee -a ${log}
    echo -e "-----Done Step 2: Lifting over from ${lift} to desired build-----\n\n" 2>&1 | tee -a ${log}
else
    echo "Assuming the file is already on desired build. No liftOver." 2>&1 | tee -a ${log}
    echo -e "-----Skipping Step 2: Lifting over from ${lift} to desired build-----\n\n" 2>&1 | tee -a ${log}
    bed=${newinfile}.bed
fi

# unifying filenames for the next step to fetch
mv ${bed} ${trait}_${BUILD}_allchr.bed

# keeping only autosomal variants and count
for i in {1..22}
do
    awk '{if ($1=="chr""'$i'") print}' ${trait}_${BUILD}_allchr.bed >> ${trait}_${BUILD}.bed
done

nauto=$(wc -l ${trait}_${BUILD}.bed | awk '{print $1}')
nall=$(wc -l ${trait}_${BUILD}_allchr.bed | awk '{print $1}')
nnauto=$((nall - nauto))
echo "Number of non-autosomal variants being discarded: $nnauto" 2>&1 | tee -a ${log}

#############################################
#@@@      Step II Running PRS per chr     @@@#
#############################################

echo "-----Step 3: Preparing to calculate PRS per chr-----" 2>&1 | tee -a ${log}
echo "Will find overlapped loci in pvar files, remove at/cg, allele-mismatching, and flip strand accordingly. " 2>&1 | tee -a ${log}

touch score.list
echo -e 'CHR('${BUILD}')\tBP('${BUILD}')\tOriginalSNPID\tUpdatedSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > ${trait}_${BUILD}_noAtCg_flipped.list
echo -e 'CHR('${BUILD}')\tBP('${BUILD}')\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > ${trait}_${BUILD}_noAtCg_mismatch.list
echo -e 'CHR('${BUILD}')\tBP('${BUILD}')\tOriginalSNPID\tOriginalRiskAllele\tOriginalRefAllele\tWeight' > ${trait}_${BUILD}_noAtCg_missing_in_pvar.list
echo -e 'CHR('${BUILD}')\tBP('${BUILD}')\tOriginalSNPID\tUpdatedSNPID\tUpdatedRiskAllele\tUpdatedRefAllele\tWeight' > ${trait}_${BUILD}_noAtCg_cleaned_forRecord.list

n_atcg_all=0 # keep the counts of at/cg loci
for CHR in {1..22}
do
    # subset to this chromosome
    awk 'BEGIN {OFS="\t"}{if ($1=="chr""'${CHR}'") print $0}' ${trait}_${BUILD}.bed > chr${CHR}_${trait}_${BUILD}.bed
    # skip if no variants present
    nline=$(wc -l chr${CHR}_${trait}_${BUILD}.bed | awk '{print $1}')
    if [[ ${nline} == "0" ]]
    then
        continue
    fi 
    # cp the plink file, make into pgen for the chr
    ${BINPATH}/plink2_mar --bfile ${BPATH}/${BFILE} --chr ${CHR} --make-pgen --out chr${CHR}_${BFILE}
    
    # remove ambiguous loci
    python2 ${SPATH}/atcg_bed.py chr${CHR}_${trait}_${BUILD}.bed chr${CHR}_${trait}_${BUILD}_noAtCg.bed # also output another bed chr${CHR}_${trait}_${BUILD}_at_cg.list
    bed="chr${CHR}_${trait}_${BUILD}_noAtCg.bed"
    nout=$(wc -l ${bed} | awk '{print $1}')
    n_atcg=$((nline - nout))
    n_atcg_all=$((n_atcg_all + n_atcg))
    
    
    # match against pvar file, flip strand, remove unmatched / tri-allelic codes etc., and make into a weight file for PLINK
    python2 ${SPATH}/bed2weightchr.py F ${bed} chr${CHR}_${BFILE}.pvar chr${CHR}_${trait}_${BUILD}_noAtCg_clean.weights 
    ## also output chr${CHR}_${trait}_${BUILD}_noAtCg_flipped.list,
    ## chr${CHR}_${trait}_${BUILD}_noAtCg_mismatch.list,
    ## chr${CHR}_${trait}_${BUILD}_noAtCg_missing_in_pvar.list
    ## chr${CHR}_${trait}_${BUILD}_noAtCg_cleaned_forRecord.list

    # remove duplicated SNPs with the same risk allele but differerent weights (could be an error from the original weight or liftOver)
    cut -f 1,2 chr${CHR}_${trait}_${BUILD}_noAtCg_clean.weights | sort -k 1 | uniq -d | cut -f 1 > chr${CHR}_duplicated.snps
    grep -vf chr${CHR}_duplicated.snps chr${CHR}_${trait}_${BUILD}_noAtCg_clean.weights > chr${CHR}_${trait}_${BUILD}_noAtCg_clean_noDup.weights
    weight="chr${CHR}_${trait}_${BUILD}_noAtCg_clean_noDup.weights"
     
    # again, skip if no variants present
    nline=$(wc -l ${weight} | awk '{print $1}')
    if [[ ${nline} == "0" ]]
    then
        rm chr${CHR}_${BFILE}.*
        continue
    fi 
    
    # Calculated PRS using PLINK2
    ${BINPATH}/plink2_mar --pfile chr${CHR}_${BFILE} \
                 --score ${weight} list-variants \
                 --out chr${CHR}_${trait}_prs
    
    # double checking if all weight SNPs were utilized
    nsnp_w=$(wc -l ${weight} | cut -d " " -f 1)
    nsnp_used=$(wc -l chr${CHR}_${trait}_prs.sscore.vars | cut -d " " -f 1)
    if [[ ${nsnp_w} != ${nsnp_used} ]]
    then
        echo "There still are weight SNPs being discarded in PLINK --score after being cleaned on CHR${CHR}. Can check the log file for trouble shooting." 2>&1 | tee -a ${log}
        #exit 1
    fi
    
    # remove the pgen file for step 3
    rm chr${CHR}_${BFILE}.*
    
    # write score name to a list
    echo "chr${CHR}_${trait}_prs.sscore" >> score.list
    # write snp info (flipped, discarded) to an integrated list
    cat chr${CHR}_${trait}_${BUILD}_noAtCg_flipped.list >> ${trait}_${BUILD}_noAtCg_flipped.list
    cat chr${CHR}_${trait}_${BUILD}_noAtCg_mismatch.list >> ${trait}_${BUILD}_noAtCg_mismatch.list
    cat chr${CHR}_${trait}_${BUILD}_noAtCg_missing_in_pvar.list >> ${trait}_${BUILD}_noAtCg_missing_in_pvar.list
    cat chr${CHR}_duplicated.snps >> ${trait}_duplicated.list
    cat chr${CHR}_${trait}_${BUILD}_noAtCg_cleaned_forRecord.list >> ${trait}_${BUILD}_noAtCg_cleaned_forRecord.list
    

done
echo -e "-----Done Step 3: Preparing to calculate PRS per chr-----\n\n" 2>&1 | tee -a ${log}

###############################
#@@    Step III Integrating   @@#
###############################

echo -e "-----Step 4: Integrating scores from chromosomes.-----" 2>&1 | tee -a ${log}
nscore=$(wc -l score.list | awk '{print $1}')
if [[ ${nscore} == "0" ]]
then
    echo "No variants are left for PRS score."
    cp ${log} ${dest}/
    exit 1
fi


Rscript ${SPATH}/integrate.R score.list ${trait}_prs.sscore

# organize the log report on # of snps
nflip=$(wc -l ${trait}_${BUILD}_noAtCg_flipped.list | awk '{print $1}')
nmis=$(wc -l ${trait}_${BUILD}_noAtCg_mismatch.list | awk '{print $1}')
nmissing=$(wc -l ${trait}_${BUILD}_noAtCg_missing_in_pvar.list | awk '{print $1}')
ndup=$(wc -l ${trait}_duplicated.list | awk '{print $1}')
nused=$(wc -l ${trait}_${BUILD}_noAtCg_cleaned_forRecord.list | awk '{print $1}')
echo "Number of ambigous A/T, C/G loci that are removed: ${n_atcg_all}" | tee -a ${log}
echo "Number of weight SNPs that are flipped to the other strand: $((nflip-1))" | tee -a ${log}
echo "Number of SNPs with mismatched allele codes against pgen that are removed: $((nmis-1))" | tee -a ${log}
echo "Number of SNPs in weight file not found in pgen that are removed: $((nmissing-1))" | tee -a ${log}
echo "Number of duplicated entries with same position and alleles (and are removed): ${ndup}" | tee -a ${log}
echo "Number of final variants used for PRS: $((nused-1))" | tee -a ${log}

echo -e "-----Done Step 4: Integrating scores from chromosomes.-----\n\n" 2>&1 | tee -a ${log}
echo -e "Copying files to destinations..." 2>&1 | tee -a ${log}
# clean up 

rm score.list
rm ${infile}
rm chr*_${trait}*
rm chr*_duplicated.snps
rm ${trait}_${BUILD}_allchr.bed
rm ${trait}_${BUILD}.bed
rm ${newinfile}
