library(stringr)
library(rtracklayer)
library(tidyr)
library(dplyr)
library(ggplot2)

workDir="/Volumes/meister.data/jsemple/ChIPseq_H3K9me2_20260405"
dir.create(paste0(workDir,"/exploreChIP"))
setwd(workDir)

## all files
files<-list.files(path=paste0(workDir,"/bigwigLog2fc"), pattern="IPlog2fcInput\\.bigwig",
           full.names=T) |> str_replace("Volumes","mnt")

samples<-basename(files) |> str_split_i(pattern="\\.",1) |> str_replace("_IPlog2fcInput","")

df<-data.frame(sampleName=samples,fileName=files)
df
write.table(df,paste0(workDir,"/exploreChIP/bwListFile_all.csv"),quote=F,row.names=F,
            col.names=F,sep=",")

## average of replicates
files<-list.files(path=paste0(workDir,"/bigwigLog2fc"), pattern="_avr.bigwig",
                  full.names=T) |> str_replace("Volumes","mnt")

samples<-basename(files) |> str_split_i(pattern="\\.",1) |> str_replace("_IPlog2fcInput_avr","")

df<-data.frame(sampleName=samples,fileName=files)
df
write.table(df,paste0(workDir,"/exploreChIP/bwListFile_avrlfc.csv"),quote=F,row.names=F,
            col.names=F,sep=",")


## consensus peaks
consensusDir=paste0(workDir,"/bwa/merged_library/macs3/broad_peak/consensus/H3K9me2")
outDir<-paste0(workDir,"/exploreChIP")
dir.create(outDir, showWarnings = FALSE)

consensus<-read.delim(paste0(consensusDir,"/H3K9me2.consensus_peaks.boolean.annotatePeaks.txt"))

# remove columns with no data
consensus$Focus.Ratio.Region.Size<-NULL
consensus$Nearest.Refseq<-NULL
consensus$Nearest.Ensembl<-NULL
consensus$Gene.Alias<-NULL
consensus$Gene.Description<-NULL

head(consensus)

fcCols<-grep("\\.*fc",colnames(consensus))
consensus$allpeakAvr<-consensus[,fcCols] |>
  mutate(row = row_number()) |>
  pivot_longer(-row) |>
  separate_rows(value, sep = ";") |>
  mutate(value = as.numeric(value)) |>
  group_by(row) |>
  summarise(avg = mean(value,na.rm=T)) |> pull("avg")


forBed<-consensus[,c("chr","start","end","interval_id","allpeakAvr")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")
hist(forBed$score,breaks=50,main=paste0("Avr all peak fc  (n=",nrow(consensus),")"))
width=forBed$chromEnd-forBed$chromStart
hist(width,breaks=50,main=paste0("Peak width (n=",nrow(consensus),")"))
export(forBed,paste0(outDir,"/H3K9me2_all.consensus_peaks.bed"))




# N2 peaks
N2boolCols<-grep("N2.*bool",colnames(consensus))
inAllN2<-rowSums(consensus[,N2boolCols])==3
result<-consensus[inAllN2, ]


N2fcCols<-grep("N2_IP_REP.*fc",colnames(result))
result$N2peakAvr<-result[,N2fcCols] |>
  mutate(row = row_number()) |>
  pivot_longer(-row) |>
  separate_rows(value, sep = ";") |>
  mutate(value = as.numeric(value)) |>
  group_by(row) |>
  summarise(avg = mean(value,na.rm=T)) |> pull("avg")

numN2peaks<-colSums(result[,N2boolCols])[1]


forBed<-result[,c("chr","start","end","interval_id","N2peakAvr")]
colnames(forBed)<-c("chrom","chromStart","chromEnd","name","score")
hist(forBed$score,breaks=50,main=paste0("Avr N2 peak fc  (n=",numN2peaks,")"))
width=forBed$chromEnd-forBed$chromStart
hist(width,breaks=50,main=paste0("Peak width (n=",numN2peaks,")"))
ggplot(data.frame(width=width),aes(x=width)) +
  geom_histogram() + theme_bw()
ggplot(data.frame(width=width),aes(x=width)) +
  geom_histogram() + theme_bw() + xlim(c(0,10000))
export(forBed,paste0(outDir,"/H3K9me2_N2.consensus_peaks.bed"))
