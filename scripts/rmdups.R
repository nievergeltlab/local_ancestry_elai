args <- commandArgs(trailingOnly = TRUE)
 infile <- args[1]

 
 
  for (chr in 13 ) # c(1:22))
  {
   a = read.table(paste("tf/", infile,"_v5_",chr,".bim",sep = ""),header=FALSE,stringsAsFactors=F)
   b = a[duplicated(a$V4),]
   write.table(b,file=paste("tf/", infile,"_v5_",chr,".snpexclude",sep = ""),sep=" ",quote = FALSE,col.names = FALSE, row.names=FALSE)
   print(paste(dim(b)[1], "deleted for chr", chr))
  }
 q("no")
