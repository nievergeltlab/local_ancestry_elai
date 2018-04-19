args <- commandArgs(trailingOnly = TRUE)
 dosfile <- args[1] #Dosage data
 famfile <- args[2] # .fam file
 subjs <- args[3] # Subjects to extract
 snps <- args[4] #SNPs to extract - supp[ly
 
 library(data.table) 
 # imputed_genotypes/dos_pts_"$stu"_mix_am-qc.hg19.ch.fl.chr13_054_057.out.dosage.gz.doscnt.gz /home/cnieverg/all_samples_ibd/pre_filter/"$file".fam "$stu".aam.subjects rs115539978.tophits

 #  dosfile="imputed_genotypes/dos_pts_fscd_mix_am-qc.hg19.ch.fl.chr13_054_057.out.dosage.gz.doscnt.gz"
 # famfile="/home/cnieverg/all_samples_ibd/pre_filter/pts_fscd_mix_am-qc.fam"
 #  subjs="fscd.aam.subjects"
 #  snps="rs115539978.tophits"
   dosages = fread(paste('zcat' , dosfile),data.table=F,sep=",",header=T)
   fam = fread(famfile,data.table=F)
   subs = fread(subjs,data.table=F)
   snplist = fread(snps,data.table=F)
   
   fam$FID_IID <- paste(fam[,1],fam[,2],sep='_')
   subs$FID_IID <- paste(subs[,1],subs[,2],sep='_')
   
   names(dosages) <- c("SNP","A1","A2",fam$FID_IID) 
   
   dosages_snps <- dosages[which(dosages$SNP %in% snplist[,1]),]
   dosages_subs <- dosages_snps[,c(1,2,3,which(names(dosages_snps) %in% subs$FID_IID))]
  
   write.csv(dosages_subs,file=paste(dosfile,".aam.doscnt",sep = ""),quote = FALSE,row.names=FALSE)

   
 