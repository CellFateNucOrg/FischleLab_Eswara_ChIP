library(dplyr)
library(rtracklayer)

workDir="/Volumes/meister.data/jsemple/ChIPseq_H3K9me2_20260405"

fqlist<-read.delim(paste0(workDir,"/fastqList.txt"),header=F)

fqlist$sample<-sapply(strsplit(fqlist$V1, "/"),"[[",7)

sampleSheet<-data.frame(sample="",
                        fastq_1=fqlist$V1[seq(1,nrow(fqlist),2)],
                        fastq_2=fqlist$V1[seq(2,nrow(fqlist),2)],
                        replicate="",antibody="",
                        control="",control_replicate="")

samples<-sapply(strsplit(sampleSheet$fastq_1, "/"),"[[",7)

if(all(samples==sapply(strsplit(sampleSheet$fastq_2, "/"),"[[",7))){
  sampleSheet$sample<-gsub("_rep[0-9]+","",samples)
  sampleSheet$replicate<-gsub("rep","",sapply(strsplit(samples, "_"),"[[",3))
  # antibody
  sampleSheet$antibody[grepl("_IP",sampleSheet$sample)]<-"H3K9me2"
  sampleSheet$antibody[grepl("Trip_IP",sampleSheet$sample)]<-"H3K9me_null"
  sampleSheet$antibody[grepl("IgG_IP",sampleSheet$sample)]<-"IgG"
  # control samples
  sampleSheet$control[sampleSheet$antibody!=""]<-gsub("_IP","_input",sampleSheet$sample[sampleSheet$antibody!=""])
  sampleSheet$control[grepl("IgG_IP",sampleSheet$sample)]<-"N2_input"
  # control replicate number
  sampleSheet$control_replicate[sampleSheet$antibody!=""]<- sampleSheet$replicate[sampleSheet$antibody!=""]
}
head(sampleSheet)


write.csv(sampleSheet,paste0(workDir,"/sampleSheet.csv"), row.names=F, quote=F)


# convert blacklist
# get blacklist from:
# "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/ce11-blacklist.v2.bed.gz"
blURL<-"https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/ce11-blacklist.v2.bed.gz"
download.file(blURL,destfile=basename(blURL))

bl<-import(paste0(workDir,"/ce11-blacklist.v2.bed"), format="BED")
seqlevels(bl)<-gsub("chr","",seqlevels(bl))
export(bl, paste0(workDir,"/WBcel235-blacklist.v2.bed"))
