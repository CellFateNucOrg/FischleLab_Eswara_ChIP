library(dplyr)

workDir="/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ChIP"

fqlist<-read.delim(paste0(workDir,"/fastqList.txt"),header=F)

fqlist$sample<-sapply(strsplit(fqlist$V1, "/"),"[[",8)

sampleSheet<-data.frame(group="",
                        fastq_1=fqlist$V1[seq(1,nrow(fqlist),2)],
                        fastq_2=fqlist$V1[seq(2,nrow(fqlist),2)],
                        replicate="",antibody="",
                        control="",control_replicate="")

samples<-sapply(strsplit(sampleSheet$fastq_1, "/"),"[[",8)

if(all(samples==sapply(strsplit(sampleSheet$fastq_2, "/"),"[[",8))){
  sampleSheet$replicate<-gsub("rep","",sapply(strsplit(samples, "_"),"[[",3))
  sampleSheet$group<-gsub("_rep[0-9]+","",samples)
  sampleSheet$antibody[grepl("_IP",sampleSheet$group)]<-"H3K9me"
  sampleSheet$control[sampleSheet$antibody=="H3K9me"]<-gsub("_IP","_input",sampleSheet$group)
  sampleSheet$control_replicate[sampleSheet$antibody=="H3K9me"]<- sampleSheet$replicate[sampleSheet$antibody=="H3K9me"]
}
head(sampleSheet)

write.csv(sampleSheet,paste0(workDir,"/sampleSheet.csv"), row.names=F, quote=F)

# convert blacklist
bl<-import(paste0(workDir,"/ce11-blacklist.v2.bed"), format="BED")
seqlevels(bl)<-gsub("chr","",seqlevels(bl))
export(bl, paste0(workDir,"/WBcel235-blacklist.v2.bed"))
