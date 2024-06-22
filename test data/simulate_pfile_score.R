## Simulates plink2 pfile in both dosages and hard call genotypes
## Simulates weight input for ESCALATOR under different situations
## Outputs expected scores to be compared against.
# Author: M. Lin 
# Date: 05/16/2024

### functions ###
bin_sum <- function(x){ # sum every two adjacent numbers into one
    if(length(x)%%2 != 0){
        stop("length of x is not even")
    }
    out <- numeric(length(x)/2)
    for(i in 1:(length(x)/2)){
        out[i] <- x[2*i-1] + x[2*i]
    }
    return(out)
}

tped2mat<- function(std_a1, std_a2, freq_a1, tped){ # convert tped w/ genocodes to dosage matrix, missing genotypes filled as external freqeuncy
    # stand_a1 or a2: vector of a1 or a2 to appear in the dosage matrix, a1 first and is which the dosage is counted on
    outmat <- matrix (NA, nrow=nrow(tped), ncol =(ncol(tped)-4)/2+3)
    outmat[,1] <- tped[,2]
    outmat[,2] <- std_a1
    outmat[,3] <- std_a2
    for(i in 1:length(std_a1)){ # iterate through each row (variant)
         genotypes <- as.character(tped[i, 5:ncol(tped)])
         tmp_geno <- as.numeric(genotypes == std_a1[i])
         tmp_geno[which(genotypes=="0")] <- freq_a1[i] # fill in missing genotypes with frequency
         dosages <- bin_sum(tmp_geno)
         outmat[i, 4:(ncol(outmat))] <- dosages
    }
    return(outmat)
}

flip_strand <- function(x){ # flip strand
    if(x == "A"){
        return("T")
    }else if(x == "T"){
        return("A")
    }else if(x == "C"){
        return("G")
    }else if(x == "G"){
        return("C")
    }else{
        stop("allele code not recognized")
    }
}

### function ends ###



# initation params
seed0 <- 12300
## number of samples
n <- 100
## number of SNPs
m <- 1000
## number of score SNPs
m_score <- 200
# allele code
allele_codes <- c("A", "C", "G", "T")

### chromosome numbers
set.seed(seed0)
chrs <- sort(round(sample(c(1:22), size=m, replace=T)))
### bp positions
bps <- c(1:m)


############## Simulation starts ##############

## dosage matrix in FORMAT=1
## readble by plink v2.0: https://www.cog-genomics.org/plink/2.0/input#import_dosage
## Dosage to pfile: plink2 --import-dosage dummy_dosage.dat format=1 noheader --fam dummy_dosage.fam --map dummy_dosage.map --make-pgen --out dummy_dosage
## Genotype file based on dosage conversion - convert to hardcall genotypes with (inevitably) a little missingness: plink2 --import-dosage dummy_dosage.dat format=1 noheader --fam dummy_dosage.fam --map dummy_dosage.map --hard-call-threshold 0.2 --make-bed --out dummy_hardcall
## convert bed file to old plink text file to be readable into R: plink1.9 --bfile dummy_hardcall --recode transpose --out dummy_hardcall

# 1. simulate dosage matrix
dosage_mat <- matrix(NA, nrow = m, ncol=n + 3)

for (i in 1:m) {
    seedi <- seed0 + i
    snpid <- paste0("SNP", i)
    set.seed(seedi)
    a1a2 <- sample(allele_codes, 2, replace = FALSE)
    a1 <- a1a2[1]
    a2 <- a1a2[2]
    set.seed(seedi)
    dosages <- runif(n, 0, 2)
    dosage_mat[i,] <- c(snpid, a1, a2, as.character(dosages))
}
write.table(dosage_mat, "dummy_dosage.dat", col.names=F, row.names=F, quote=F, sep=' ')
## 2. simulate companion fam and map files
fam <- cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), rep("0", n), rep("0", n), rep("0", n), rep("0", n)) # FID, IID, pID, mID, sex, phenotype
map <- cbind(as.character(chrs), dosage_mat[,1], rep("0", m), as.character(bps)) #chr, rsID, genetic distance, bp
write.table(fam, "dummy_dosage.fam", col.names=F, row.names=F, quote=F, sep=' ')
write.table(map, "dummy_dosage.map", col.names=F, row.names=F, quote=F, sep=' ')

## 2.1 convert to pfile for plink2
### dosage
system('plink2 --import-dosage dummy_dosage.dat format=1 noheader --fam dummy_dosage.fam --map dummy_dosage.map --make-pgen --out dummy_dosage')
### hard call genotypes (will have missingness)
system('plink2 --import-dosage dummy_dosage.dat format=1 noheader --fam dummy_dosage.fam --map dummy_dosage.map --hard-call-threshold 0.2 --make-bed --out dummy_hardcall')
system('plink2 --bfile dummy_hardcall --make-pgen --out dummy_hardcall')




## 3. simulate weight files (for different levels of testing purposes, won't cover liftOver if using simulated data)
### v0: all weights are 1, no strand flipping etc.
cmd0 <- 'echo "#genome_build = GRCh38" > dummy_lv0.wgt'
cmd1 <- 'echo -e "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" >> dummy_lv0.wgt'
set.seed(seed0)
index <- sort(sample(c(1:m), m_score, replace = FALSE)) # which SNPs to use as score SNPs
weight0 <- rep(1, m_score)
weight0_mat <- cbind(chrs[index], bps[index], dosage_mat[index,2], dosage_mat[index,3], as.character(weight0))
write.table(weight0_mat, "dummy_lv0.wgt_tmp", col.names=F, row.names=F, quote=F, sep='\t')
cmd2 <- 'cat dummy_lv0.wgt_tmp >> dummy_lv0.wgt'
cmd3 <- 'rm dummy_lv0.wgt_tmp'
system(cmd0)
system(cmd1)
system(cmd2)
system(cmd3)

### v1: all weights are random, no other issue
cmd0 <- 'echo "#genome_build = GRCh38" > dummy_lv1.wgt'
cmd1 <- 'echo -e "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" >> dummy_lv1.wgt'
set.seed(seed0)
weight1 <- runif(m_score, -5, 10)
weight1_mat <- cbind(chrs[index], bps[index], dosage_mat[index,2], dosage_mat[index,3], as.character(weight1))
write.table(weight1_mat, "dummy_lv1.wgt_tmp", col.names=F, row.names=F, quote=F, sep='\t')
cmd2 <- 'cat dummy_lv1.wgt_tmp >> dummy_lv1.wgt'
cmd3 <- 'rm dummy_lv1.wgt_tmp'
system(cmd0)
system(cmd1)
system(cmd2)
system(cmd3)

### v2. use weights from v1, but skip one chromosome (chr2)
cmd0 <- 'echo "#genome_build = GRCh38" > dummy_lv2.wgt'
cmd1 <- 'echo -e "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" >> dummy_lv2.wgt'
which_chr2 <- which(chrs[index] == 2)
index_nochr2 <- index[-which_chr2]
weight2_mat <- cbind(chrs[index_nochr2], bps[index_nochr2], dosage_mat[index_nochr2,2], dosage_mat[index_nochr2,3], as.character(weight1[-which_chr2]))
write.table(weight2_mat, "dummy_lv2.wgt_tmp", col.names=F, row.names=F, quote=F, sep='\t')
cmd2 <- 'cat dummy_lv2.wgt_tmp >> dummy_lv2.wgt'
cmd3 <- 'rm dummy_lv2.wgt_tmp'
system(cmd0)
system(cmd1)
system(cmd2)
system(cmd3)

### v3. use weights from v1, with some strand flipping: avoid including ambigous variants
cmd0 <- 'echo "#genome_build = GRCh38" > dummy_lv3.wgt'
cmd1 <- 'echo -e "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" >> dummy_lv3.wgt'
score_alleles <- paste0(dosage_mat[index,2], ":", dosage_mat[index,3])
which_amb <- which(score_alleles%in% c("A:T", "T:A", "C:G", "G:C")) # which of the lv1 variants are ambigous
set.seed(seed0)
which_flip <- sample(c(1:m_score)[-which_amb], 5, replace = FALSE) # make sure those to flip are not ambigous coded
a1_v3 <- dosage_mat[index, 2] # original a1 code
a1_v3[which_flip] <- as.character(sapply(a1_v3[which_flip], flip_strand)) # flip strand 5 SNPs
a2_v3 <- dosage_mat[index, 3] # original a2 code 
a2_v3[which_flip] <- as.character(sapply(a2_v3[which_flip], flip_strand)) # flip strand 5 SNPs
weight3_mat <- cbind(chrs[index], bps[index], a1_v3, a2_v3, as.character(weight1))
write.table(weight3_mat, "dummy_lv3.wgt_tmp", col.names=F, row.names=F, quote=F, sep='\t')
cmd2 <- 'cat dummy_lv3.wgt_tmp >> dummy_lv3.wgt'
cmd3 <- 'rm dummy_lv3.wgt_tmp'
system(cmd0)
system(cmd1)
system(cmd2)
system(cmd3)

## v4. use weights from v1, with some non-matching allele code, and extra weight SNPs not found
cmd0 <- 'echo "#genome_build = GRCh38" > dummy_lv4.wgt'
cmd1 <- 'echo -e "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" >> dummy_lv4.wgt'
### add SNPs in the weight with non-matching allele codes
set.seed(seed0)
index_nonmatch <- sort(sample(c(1:m)[-index], 10, replace = FALSE)) # pick 10 SNPs that are not in the causal list, to add as additional var in the weight as non-matching allele codes
a1_all <- dosage_mat[, 2] # all original a1 code
a2_all <- dosage_mat[, 3] # all original a2 code
count <- 0 # count how many SNPs have an allele code been replaced, total 10
for(j in index_nonmatch){
    count <- count + 1
    a1 <- a1_all[j]
    a2 <- a2_all[j]
    third_allele <- setdiff(c("A", "C", "G", "T"), c(a1, a2))[1]
    if(count <=5){ # first 5 SNPs, replace a1
        a1_all[j] <- third_allele
    }else{ # last 5 SNPs, replace a2
        a2_all[j] <- third_allele
    }
}
set.seed(seed0)
weight_nonmatch <- runif(10, min=11, max=20)
index_v4_sub <- sort(c(index, index_nonmatch)) # index of SNPs in the v1 and those with non-matching allele codes
temp_reposition_index <- match(index_v4_sub, c(index, index_nonmatch))
weight4_sub <- c(weight1, weight_nonmatch)[temp_reposition_index] # weight for all SNPs in the v1 and those with non-matching allele codes, in the same order as in their sorted chr:bp

weight4_mat_sub <- cbind(chrs[index_v4_sub], bps[index_v4_sub], a1_all[index_v4_sub], a2_all[index_v4_sub], as.character(weight4_sub))
### add SNPs in the weight with positions not found in the simulated data
set.seed(seed0)
chr_extra <- sample(c(1:22), 5, replace = FALSE)
set.seed(seed0)
bp_extra <- sample(setdiff(c(1:10000), bps), 5, replace = FALSE) # 5 random positions not in the simulated data
set.seed(seed0)
weight_extra <- runif(5, min=21, max=30)
set.seed(seed0)
weight4_mat_extra <- cbind(chr_extra, bp_extra, rep("A", 5), rep("G", 5), as.character(weight_extra))

weight4_mat_unsorted <- rbind(weight4_mat_sub, weight4_mat_extra)
temp_v4_chr_bp <- as.numeric(weight4_mat_unsorted[,1])*1e9 + as.numeric(weight4_mat_unsorted[,2]) # combine chr and bp to sort
index_v4_mat <- sort(temp_v4_chr_bp, index.return=T)$ix
weight4_mat <- weight4_mat_unsorted[index_v4_mat,]

write.table(weight4_mat, "dummy_lv4.wgt_tmp", col.names=F, row.names=F, quote=F, sep='\t')
cmd2 <- 'cat dummy_lv4.wgt_tmp >> dummy_lv4.wgt'
cmd3 <- 'rm dummy_lv4.wgt_tmp'
system(cmd0)
system(cmd1)
system(cmd2)
system(cmd3)

### v5. (complicated situation) v2-v4 combined
#### use index from v4 since it's the largest
#### add in flipped index and remove chr2 (avoid ambigous variants to flip)
cmd0 <- 'echo "#genome_build = GRCh38" > dummy_lv5.wgt'
cmd1 <- 'echo -e "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" >> dummy_lv5.wgt'
score_alleles_v4 <- paste0(weight4_mat[,3], ":", weight4_mat[,4])
which_amb_v4 <- which(score_alleles_v4%in% c("A:T", "T:A", "C:G", "G:C")) # which of the lv4 variants are ambigous
set.seed(seed0)
which_v5_flip <- sample(c(1:dim(weight4_mat)[1])[-which_amb_v4], size=5, replace = FALSE) # only flip those that are not ambigous
a1_v5 <- weight4_mat[,3] 
a2_v5 <- weight4_mat[,4]
a1_v5[which_v5_flip] <- as.character(sapply(a1_v5[which_v5_flip], flip_strand)) # flip strand 5 SNPs
a2_v5[which_v5_flip] <- as.character(sapply(a2_v5[which_v5_flip], flip_strand)) # flip strand 5 SNPs
weight5_mat_temp <- weight4_mat
weight5_mat_temp[,3] <- a1_v5
weight5_mat_temp[,4] <- a2_v5
which_v5_chr2 <- which(weight5_mat_temp[,1] == 2)
weight5_mat <- weight5_mat_temp[-which_v5_chr2,]
write.table(weight5_mat, "dummy_lv5.wgt_tmp", col.names=F, row.names=F, quote=F, sep='\t')
cmd2 <- 'cat dummy_lv5.wgt_tmp >> dummy_lv5.wgt'
cmd3 <- 'rm dummy_lv5.wgt_tmp'
system(cmd0)
system(cmd1)
system(cmd2)
system(cmd3)



## 4. Generate the "correct" scores to compare with those from ESCALATOR

### 4.1 first read in the hardcall genotypes in tped and convert, along with freq file to fill in imputation
system('plink1.9 --bfile dummy_hardcall --recode transpose --out dummy_hardcall')
tped <- read.table("dummy_hardcall.tped")
#### also need frequency file from PLINK to fill in missing genotypes: 
system('plink2 --bfile dummy_hardcall --freq --out dummy_hardcall')
freqfile <- read.table("dummy_hardcall.afreq")
freqs <- freqfile[,5]
index_change_freq <- which(freqfile[,4] != dosage_mat[,2]) # mark those where alt allele in freq file is not the a1 file in dosage mat
freqs[index_change_freq] <- 1 - freqs[index_change_freq] # flip the frequency
#### convert to dosage matrix
tped_mat <- tped2mat(dosage_mat[,2], dosage_mat[,3], freqs, tped) 
tped_mat_sub <- apply(tped_mat[index, c(4:(n+3))],2, as.numeric)# subset to score SNPs, this will be used to calculate scores

### 4.2 genenate score scores when ambigous variants are NOT removed
dosage_mat_sub <- apply(dosage_mat[index,c(4:(n+3))],2, as.numeric)
#### lv0
correct_score_lv0_dosage <- weight0%*%dosage_mat_sub # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv0_dosage)), "dummy_lv0_correct_wAmb_dosage_based.score", col.names=F, row.names=F, quote=F, sep=' ')
correct_score_lv0_geno <- weight0%*%tped_mat_sub # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv0_geno)), "dummy_lv0_correct_wAmb_geno_based.score", col.names=F, row.names=F, quote=F, sep=' ')
#### lv1 = lv3 = lv4
correct_score_lv1_dosage <- weight1%*%dosage_mat_sub # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv1_dosage)), "dummy_lv1_correct_wAmb_dosage_based.score", col.names=F, row.names=F, quote=F, sep=' ')
correct_score_lv1_geno <- weight1%*%tped_mat_sub # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv1_geno)), "dummy_lv1_correct_wAmb_geno_based.score", col.names=F, row.names=F, quote=F, sep=' ')
#### lv2 = lv5
weight2 <- weight1[-which_chr2]
dosage_mat_sub_v25 <- dosage_mat_sub[-which_chr2,]
tped_mat_sub_v25 <- tped_mat_sub[-which_chr2,]
correct_score_lv2_dosage <- weight2%*%dosage_mat_sub_v25 # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv2_dosage)), "dummy_lv2_correct_wAmb_dosage_based.score", col.names=F, row.names=F, quote=F, sep=' ')
correct_score_lv2_geno <- weight2%*%tped_mat_sub_v25 # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv2_geno)), "dummy_lv2_correct_wAmb_geno_based.score", col.names=F, row.names=F, quote=F, sep=' ')


### 4.3 generate scores when ambigous variants are removed
#score_alleles <- paste0(dosage_mat[index,2], ":", dosage_mat[index,3])
#which_amb <- which(score_alleles%in% c("A:T", "T:A", "C:G", "G:C"))
index_nonamb <- index[-which_amb] # index in the entire m variant of non-ambigous variants
weight0_nonamb <- weight0[-which_amb]
weight1_nonamb <- weight1[-which_amb]
which_chr2_nonamb <- which(chrs[index_nonamb] == 2)
weight2_nonamb <- weight1_nonamb[-which_chr2_nonamb]
dosage_mat_sub_nonamb <- apply(dosage_mat[index_nonamb,c(4:(n+3))],2, as.numeric)
tped_mat_sub_nonamb <- apply(tped_mat[index_nonamb, c(4:(n+3))], 2, as.numeric) # subset to score SNPs, this will be used to calculate scores
#### lv0
correct_score_lv0_dosage_nonamb <- weight0_nonamb%*%dosage_mat_sub_nonamb # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv0_dosage_nonamb)), "dummy_lv0_correct_woAmb_dosage_based.score", col.names=F, row.names=F, quote=F, sep=' ')
correct_score_lv0_geno_nonamb <- weight0_nonamb%*%tped_mat_sub_nonamb # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv0_geno_nonamb)), "dummy_lv0_correct_woAmb_geno_based.score", col.names=F, row.names=F, quote=F, sep=' ')
#### lv1 = lv3 = lv4
correct_score_lv1_dosage_nonamb <- weight1_nonamb%*%dosage_mat_sub_nonamb # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv1_dosage_nonamb)), "dummy_lv1_correct_woAmb_dosage_based.score", col.names=F, row.names=F, quote=F, sep=' ')
correct_score_lv1_geno_nonamb <- weight1_nonamb%*%tped_mat_sub_nonamb # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv1_geno_nonamb)), "dummy_lv1_correct_woAmb_geno_based.score", col.names=F, row.names=F, quote=F, sep=' ')
#### lv2 = lv5
dosage_mat_sub_nonamb_v25 <- dosage_mat_sub_nonamb[-which_chr2_nonamb,]
tped_mat_sub_nonamb_v25 <- tped_mat_sub_nonamb[-which_chr2_nonamb,]
correct_score_lv2_dosage_nonamb <- weight2_nonamb%*%dosage_mat_sub_nonamb_v25 # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv2_dosage_nonamb)), "dummy_lv2_correct_woAmb_dosage_based.score", col.names=F, row.names=F, quote=F, sep=' ')
correct_score_lv2_geno_nonamb <- weight2_nonamb%*%tped_mat_sub_nonamb_v25 # correct score to compare
write.table(cbind(paste0("F", c(1:n)), paste0("I", c(1:n)), as.character(correct_score_lv2_geno_nonamb)), "dummy_lv2_correct_woAmb_geno_based.score", col.names=F, row.names=F, quote=F, sep=' ')



### 5. parse chr and clean up, and store the object
system('for i in {1..22}; do plink2 --pfile dummy_dosage --chr $i --make-pgen --out chr${i}_dummy_dosage; done')
system('for i in {1..22}; do plink2 --pfile dummy_hardcall --chr $i --make-pgen --out chr${i}_dummy_hardcall; done')
#system('rm dummy_dosage.*')
#system('rm dummy_hardcall.*')
save.image('simulate_dummy_data.RData')
