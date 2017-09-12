cat results/chr*_caseonly_samplesite.results | grep -v regsums > results/all_chr.results

library(data.table)
results <- fread('results/all_chr.results', header=F,data.table=F)
names(results) <- c("rs","pos","chr","maf","P")

#Get Bayes values    
posterior <- function(x,prior,lambda) {(dchisq(x,1,lambda)*prior)/((dchisq(x,1,lambda)*prior)+(dchisq(x,1,0)*(1-prior)))}
admixture_burden <- 202 #Estimated 6/27, may need to fix some chrs, estimated using the 03_estimate_prior script

admixture_lambda <- (qnorm(1-0.05/admixture_burden/2)+qnorm(0.8))^2
admixture_prior <- 1/admixture_burden
admixture_p <- results$P

admixture_test <- qchisq(admixture_p,1,0,lower.tail=FALSE)
admixture_posterior <- posterior(x=admixture_test,prior=admixture_prior,lambda=admixture_lambda)


results$bf <- admixture_posterior


source('manhattan_bayes.R')
source('manhattan_plot.R')

ManhattanPlot_bayes(results$chr, results$pos, results$bf, results$rs, genomesig = .5, genomesug = 0,photoname = '', 
    outname = 'results/admixture_mh_plot.jpg', colors = c(rgb(230,185,184,max=255),'black'), sigcolor = 'darkred', sugcolor = 'indianred', ncex = .5, 
    plotsymbol = 16, sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA',ylabel = "Posterior Probability of Association")

    
    ManhattanPlot_AJS_cut(results$chr, results$pos, results$P, results$rs, genomesig = .00025, genomesug = 0,photoname = '', 
    outname = 'results/admixture_mh_plot_traditional.jpg', colors = c(rgb(230,185,184,max=255),'black'), sigcolor = 'darkred', sugcolor = 'indianred', ncex = .5, 
    plotsymbol = 16, sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA')

    
#Annotation of singificant data
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # for annotation
library(org.Hs.eg.db)     


sigdata <- subset(results, bf > .5,select=c(rs,chr,pos))
names(sigdata) <- c("rsid","chr","pos")

sigdata$chr <- paste("chr",sigdata$chr,sep='')

target <- with(sigdata,
                GRanges( seqnames = Rle(chr),
                         ranges   = IRanges(pos, end=pos, names=rsid),
                         strand   = Rle(strand("*")) ) )

                         
loc <- locateVariants(target, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
names(loc) <- NULL
out <- as.data.frame(loc)
out$names <- names(target)[ out$QUERYID ]
out <- out[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")]
out <- unique(out)
out

Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)
 
x <- unique( with(out, c(levels(GENEID), levels(PRECEDEID), levels(FOLLOWID))) )
table( x %in% names(id2Symbol) ) # good, all found
 
out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
out$PRECEDESYMBOL <- id2Symbol[ as.character(out$PRECEDEID) ]
out$FOLLOWSYMBOL <- id2Symbol[ as.character(out$FOLLOWID) ]
out

out2 <- out
names(out2)[1] <- "rs"
out3 <- merge(out2,results,by="rs")


write.table(names(table(subset(out3,bf > .5)$GENESYMBOL)), 'glaucoma_genes.txt')

  
#Compare to Khor 2016 PACG
khor <- read.table('Khor-27064256.txt',header=T,stringsAsFactors=F,nr=800000)
names(khor)[1] <- 'rs'

compare <- merge(out3,khor,by="rs",suffixes=c("_ours","khors"))

compare_replicate <- subset(compare,bf>=.5)

table(compare_replicate[which(compare_replicate$Pkhors < .0005),]$GENESYMBOL )

#Compare to Bailey 2016 POAG

bailey <- fread('SUMMARY_STATS_NG_MANUSCRIPT_allchr.txt', header=F,data.table=F)

names(bailey) <- c('CHR','BP','rs','P')

compare <- merge(out3,bailey,by="rs",suffixes=c("_ours","bailey"))

compare_replicate <- subset(compare,bf>=.5)

table(compare_replicate[which(compare_replicate$Pbailey < .05),]$GENESYMBOL )


