#! /usr/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=4


source $CONDA_ACTIVATE deeptools

workDir=/mnt/external.data/MeisterLab/FischleLab_KarthikEswara/ChIP
samplesheet=$workDir/sampleSheet.csv
#sizeFactorFile=$workDir/$runName/other/deseq2/HPL2GFP_lin61_vs_N2.deseq2.sizefactors.tsv
blackListFile=$workDir/WBcel235-blacklist.v2.bed
outDirAvr=$workDir/bigwigAvrSeparate
outDirSubtract=$workDir/bigwigSubtract
outDirLog2fc=$workDir/bigwigLog2fc

bigwigPath=$workDir/bwa/merged_library/bigwig

mkdir -p $outDirAvr
mkdir -p $outDirSubtract
mkdir -p $outDirLog2fc


#### bigwigAverage ######


groups=(`tail -n +2 "$samplesheet" | cut -f1 -d"," | grep "_IP" | uniq`)
samples=("${groups[@]/%_IP/}")

for g in "${!groups[@]}"; do
  group=${groups[g]}
  sample=${samples[g]}
  echo "ips in group ${group}:"
  ips=(`awk -F',' -v val="$group" 'NR>1 && $1 == val { print $1 "_REP" $4 }' "$samplesheet"`)
  echo "${ips[@]}"
  inputs=(`awk -F',' -v val="$group" 'NR>1 && $1 == val { print $6 "_REP" $4 }' "$samplesheet"`)
  
  # add filename suffix
  ipfiles=("${ips[@]/%/.mLb.clN.bigWig}")
  inputfiles=("${inputs[@]/%/.mLb.clN.bigWig}")
  
  # add path prefix
  ipfiles=("${ipfiles[@]/#/${bigwigPath}\/}")
  inputfiles=("${inputfiles[@]/#/${bigwigPath}\/}")
  
  echo "averaging ips"
  bigwigAverage -b ${ipfiles[@]} -bs 10 -bl $blackListFile -p $SLURM_NTASKS -o ${outDirAvr}/${sample}_IP_avr.bigwig
  
  echo "averaging inputs"
  bigwigAverage -b ${inputfiles[@]} -bs 10 -bl $blackListFile -p $SLURM_NTASKS -o ${outDirAvr}/${sample}_Input_avr.bigwig

 for i in "${!ips[@]}"; do
	echo "subtracting input from ip"
	bigwigCompare --bigwig1 ${ipfiles[i]} --bigwig2 ${inputfiles[i]} \
	--operation subtract --binSize 10 -p $SLURM_NTASKS -bl $blackListFile \
	-o ${outDirSubtract}/${ips[i]}_IPminusInput.bigwig

	echo "log2 ratio of ip/input"
	bigwigCompare --bigwig1 ${ipfiles[i]} --bigwig2 ${inputfiles[i]} \
	--operation log2 --binSize 10 -p $SLURM_NTASKS -bl $blackListFile -o ${outDirLog2fc}/${ips[i]}_IPlog2fcInput.bigwig
 done
  
  echo "averaging ip minus input tracks"
  # add filename suffix
  ipfiles=("${ips[@]/%/_IPminusInput.bigwig}")
  # add path prefix
  ipfiles=("${ipfiles[@]/#/${outDirSubtract}\/}")

  bigwigAverage -b ${ipfiles[@]} -bs 10 -bl $blackListFile -p $SLURM_NTASKS  -o ${outDirSubtract}/${sample}_IPminusInput_avr.bigwig
  
  echo "averaging log2(ip/input) tracks"
  # add filename suffix
  ipfiles=("${ips[@]/%/_IPlog2fcInput.bigwig}")
  # add path prefix
  ipfiles=("${ipfiles[@]/#/${outDirLog2fc}\/}")

  bigwigAverage -b ${ipfiles[@]} -bs 10 -bl $blackListFile -p $SLURM_NTASKS  -o ${outDirLog2fc}/${sample}_IPlog2fcInput_avr.bigwig
done