args = commandArgs(T)
score.list <- read.table(args[1])
output.name <- args[2]


i=0
for(filename in as.character(score.list$V1)){
  #infile <- read.table(paste("chr", i, "_CCPM_trait_snps.sscore", sep=""), colClasses = c("character", "numeric", "numeric", "numeric"))
  i = i+1
  infile <- read.table(filename, colClasses = c("character", "numeric", "numeric", "numeric"))
  if(i==1){
    scores <- data.frame(ID=as.character(infile$V1), scores=as.numeric(infile$V4))
  }else{
    temp_score <- infile$V4[match(scores$ID, infile$V1)]
    temp_score2 <- temp_score + scores$scores
    scores$scores <- as.numeric(temp_score2)
  }
}

write.table(scores, output.name, col.names=F, row.names=F, quote=F, sep='\t')
