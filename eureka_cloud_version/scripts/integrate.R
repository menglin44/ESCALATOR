args = commandArgs(T)
score.list <- read.table(args[1])
output.name <- args[2]


i=0
for(filename in as.character(score.list$V1)){
  #infile <- read.table(paste("chr", i, "_CCPM_trait_snps.sscore", sep=""), colClasses = c("character", "numeric", "numeric", "numeric"))
  i = i+1
  #infile <- read.table(filename, colClasses = c("character", "numeric", "numeric", "numeric"))
  infile <- read.table(filename, colClasses="character")
  coln <- ncol(infile)
  if(i==1){
    if(coln==5){
      # has FID present
      scores <- data.frame(FID=as.character(infile$V1), IID=as.character(infile$V2), scores=as.numeric(infile$V5))
      scores$tempID <- paste(scores$FID, scores$IID, sep=':')
    }else if(coln==4){
      #only ID present 
    scores <- data.frame(ID=as.character(infile$V1), scores=as.numeric(infile$V4))
    }else{
      stop(paste0('Wrong number of columns in the score output for ', filename))
      }
  }else{
    if(coln==5){
      infile$tempID <- paste(infile$V1, infile$V2, sep=':')
      temp_score <- infile$V5[match(scores$tempID, infile$tempID)]
      temp_score2 <- temp_score + scores$scores
      scores$scores <- as.numeric(temp_score2)
    }else if(coln==4){
      temp_score <- infile$V4[match(scores$ID, infile$V1)]
      temp_score2 <- temp_score + scores$scores
      scores$scores <- as.numeric(temp_score2)
    }else{
      stop(paste0('Wrong number of columns in the score output for ', filename))
      }
  }
}

if(ncol(scores)==3){
  scores2 <- scores[,-1] # trim the last column of tempID
}else{
  scores2 <- scores
}

write.table(scores2, output.name, col.names=F, row.names=F, quote=F, sep='\t')
