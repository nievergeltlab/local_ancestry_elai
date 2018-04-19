args <- commandArgs(trailingOnly = TRUE)
 infile <- args[1]

 
   dat1 = read.table(infile,header=FALSE,stringsAsFactors=F)
   names(dat1) <- c("FID","IID","FID2","IID2")
   dat1$FID2_rename <- dat1$FID2 #Make a new entry of what the FID should be called

   try(dat1[duplicated(dat1$FID2_rename),]$FID2_rename <- paste(dat1[duplicated(dat1$FID2_rename),]$FID2_rename,"_duplicate",sep=""))
   try(dat1[duplicated(dat1$FID2_rename),]$FID2_rename <- paste(dat1[duplicated(dat1$FID2_rename),]$FID2_rename,"_duplicate",sep=""))
   try(dat1[duplicated(dat1$FID2_rename),]$FID2_rename <- paste(dat1[duplicated(dat1$FID2_rename),]$FID2_rename,"_duplicate",sep=""))

   dat1$IID2_rename <- dat1$FID2_rename
   print("Writing updated IDs to file")
   print(paste(infile,"_rn",sep = ""))
   write.table(subset(dat1,select=c(FID,IID,FID2_rename,IID2_rename)),file=paste(infile,"_rn",sep = ""),quote = FALSE,col.names = FALSE, row.names=FALSE)

 