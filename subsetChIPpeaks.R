library(rtracklayer)
library(GenomicRanges)


workDir="/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ChIP"
setwd(workDir)

consensusDir=paste0(workDir,"/bwa/merged_library/macs3/broad_peak/consensus/H3K9me")
outDir<-paste0(workDir,"/consensusPeaks_subsets")
dir.create(outDir, showWarnings = FALSE)

consensus<-read.delim(paste0(consensusDir,"/H3K9me.consensus_peaks.boolean.annotatePeaks.txt"))

consensus$Focus.Ratio.Region.Size<-NULL
consensus$Nearest.Refseq<-NULL
consensus$Nearest.Ensembl<-NULL
consensus$Gene.Alias<-NULL
consensus$Gene.Description<-NULL

head(consensus)
boolCols<-colnames(consensus)[grep(".bool",colnames(consensus))]

# all but Trip
boolPattern<-!grepl("Trip",boolCols)

matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_allButTrip.consensus_peaks.bed"))
#write.table(result, paste0(outDir,"/H3K9me2.consensus_peaks.boolean.annotatePeaks.allButTrip.txt"), row.names=F, quote=F, sep="\t")


# EM91 only
boolPattern<-grepl("EM91",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_EM91only.consensus_peaks.bed"))


# Trip only
boolPattern<-grepl("Trip",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_Triponly.consensus_peaks.bed"))

# all
boolPattern<-rep(TRUE, length(boolCols))
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_all.consensus_peaks.bed"))

# all but EM91 and Trip
boolPattern<-!grepl("Trip",boolCols) & !grepl("EM91",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_allButTripEM91.consensus_peaks.bed"))


# N2 only
boolPattern<-grepl("N2",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_N2only.consensus_peaks.bed"))


# EM92 only
boolPattern<-grepl("EM92",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_EM92only.consensus_peaks.bed"))


# N2 and EM91 only
boolPattern<-grepl("N2",boolCols) | grepl("EM91",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_N2&EM91only.consensus_peaks.bed"))



# EM88, EM90 and EM92 only
boolPattern<-grepl("EM88",boolCols) | grepl("EM90",boolCols) | grepl("EM92",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_EM88-EM90-EM92only.consensus_peaks.bed"))



# EM88, EM90 ,EM91 and EM92 only
boolPattern<-grepl("EM88",boolCols) | grepl("EM90",boolCols) | grepl("EM91",boolCols) | grepl("EM92",boolCols)
boolPattern
matches <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
result  <-consensus[matches, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_EM88-EM90-EM91-EM92only.consensus_peaks.bed"))


### final selection
# all but Trip
boolPattern<-!grepl("Trip",boolCols)
boolPattern
matches1 <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
# all
boolPattern<-rep(TRUE, length(boolCols))
boolPattern
matches2 <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)
# all but EM91 and Trip
boolPattern<-!grepl("Trip",boolCols) & !grepl("EM91",boolCols)
boolPattern
matches3 <- rowSums(consensus[boolCols] == matrix(boolPattern, nrow(consensus), length(boolPattern), byrow=TRUE)) == length(boolPattern)

result  <-consensus[matches1 | matches2 | matches3, ]

colSums(result[,boolCols])

forBed<-result[,c("chr","start","end","interval_id","num_peaks")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")

export(forBed,paste0(outDir,"/H3K9me2_finalSubset.consensus_peaks.bed"))
write.table(result, paste0(outDir,"/H3K9me2.consensus_peaks.boolean.annotatePeaks.finalSubset.txt"), row.names=F, quote=F, sep="\t")


gr<-GRanges(consensus)
EM88sig<-import(paste0(workDir,"/bigwigLog2fc/EM88_IPlog2fcInput_avr.bigwig"),which=gr,
               as="Numeric")

gr$meanH3K9me3<-sapply(EM88sig,mean,na.rm=T)
hist(gr$meanH3K9me3, breaks=100)

grsub<-GRanges(forBed)
EM88sig<-import(paste0(workDir,"/bigwigLog2fc/EM88_IPlog2fcInput_avr.bigwig"),which=grsub,
                as="Numeric")

grsub$meanH3K9me3<-sapply(EM88sig,mean,na.rm=T)
hist(grsub$meanH3K9me3, breaks=100)

grNotSub<-gr[! (gr$interval_id %in% grsub$name)]

EM88sig<-import(paste0(workDir,"/bigwigLog2fc/EM88_IPlog2fcInput_avr.bigwig"),which=grNotSub, as="Numeric")

grNotSub$meanH3K9me3<-sapply(EM88sig,mean,na.rm=T)
hist(grNotSub$meanH3K9me3, breaks=100)
