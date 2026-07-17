#! /bin/bash
#SBATCH --time=0-08:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=4

source $CONDA_ACTIVATE deeptools

workDir=/mnt/meister.data/jsemple/ChIPseq_H3K9me2_20260405

chipData="H3K9me2_allRep"
bwlistfile=${workDir}/exploreChIP/bwListFile_all.csv

#chipData="H3K9me2_avrlfc"
#bwlistfile=${workDir}/exploreChIP/bwListFile_avrlfc.csv

outDir=${workDir}/exploreChIP
mkdir -p $outDir




sampleNames=( $(cut -f 1 -d',' $bwlistfile) )
echo ${sampleNames[@]}
fileNames=( $(cut -f 2 -d',' $bwlistfile) )
echo ${fileNames[@]}


## bigwig files
bwPaths=${fileNames[@]}
bwNames=${sampleNames[@]}

## bed files 
peaks=( ${workDir}/exploreChIP/H3K9me2_N2.consensus_peaks.bed ${workDir}/exploreChIP/H3K9me2_all.consensus_peaks.bed )
peakNames=( "N2" "All" )

peakWidth=2000
updown=2000
regionWidths=${peakWidth}_ud${updown}

echo "number of tasks: " $SLURM_NTASKS

## compute matrix
redoMatrix=true
if $redoMatrix; then 
      blacklist=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz
      blacklistFile=`basename ${blacklist%.gz}`
      if [ ! -f "$blacklistFile" ]; then
            wget $blacklist
            gunzip `basename $blacklist`
      fi
 
      
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${peaks[0]} \
            --beforeRegionStartLength 2000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 2000 \
            --skipZeros \
            -o ${outDir}/matrix_${peakNames[0]}_${chipData}_${regionWidths}.mat.gz \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose
            #--outFileNameMatrix ${outDir}/matrix_${peakNames[0]}_${chipData}.tab \
      
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${peaks[1]} \
            --beforeRegionStartLength 2000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 2000 \
            --skipZeros \
            -o ${outDir}/matrix_${peakNames[1]}_${chipData}_${regionWidths}.mat.gz \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose
	    #--outFileNameMatrix ${outDir}/matrix_${peakNames[1]}_${chipData}.tab \
fi

# multiBigwigSummary BED-file --bwfiles ${bwPaths[@]} --BED $activeEnhancers_fountain \
#     --outRawCounts matrix_enhVhistoneModEncode_avrRaw.tab -o matrix_enhVhistoneModEncode_avrRaw.npz \
#     --numberOfProcessors 4 

#maxs=( `awk 'NR > 1 { for (i = 4; i <= NF; i++) { if ($i > max[i]) max[i] = $i+0 } } END { for (i = 4; i <= NF; i++) print max[i]/2+0 }' ${outDir}/matrix_${peakNames[0]}_${chipData}.tab` )
#maxs=( 30 `printf ' 2 %.0s' {1..21}` )
#mins=( `printf ' -2 %.0s' {1..22}` )
#maxs=( 60 300 1000 500 10 200 800 200 1000 )
#mins=( `printf ' 0 %.0s' {1..9}` )


## make plots
plotHeatmap -m ${outDir}/matrix_${peakNames[0]}_${chipData}_${regionWidths}.mat.gz \
      -out ${outDir}/${peakNames[0]}_${chipData}_heatmap_${regionWidths}_N2sort.pdf \
      --colorMap Blues  \
      --startLabel "" --endLabel "" \
      -y "" -x "Distance" \
      --regionsLabel ${peakNames[0]}  \
      --plotTitle "N2 consensus peaks" \
      --samplesLabel ${bwNames[@]} \
      --sortUsingSamples 1 2 3
#      --sortRegions no \
#      --zMin ${mins[@]} \
#      --zMax ${maxs[@]} \
#      --yMin ${mins[@]} \
#      --yMax ${maxs[@]} 

## make plots
plotHeatmap -m ${outDir}/matrix_${peakNames[1]}_${chipData}_${regionWidths}.mat.gz \
      -out ${outDir}/${peakNames[1]}_${chipData}_heatmap_${regionWidths}_N2sort.pdf \
      --colorMap Blues  \
      --startLabel "" --endLabel "" \
      -y "" -x "Distance" \
      --regionsLabel ${peakNames[1]} \
      --plotTitle "All consensus peaks" \
      --samplesLabel ${bwNames[@]} \
      --sortUsingSamples 1 2 3
      #--sortRegions no \
      #--zMin ${mins[@]} \
      #--zMax ${maxs[@]} \
      #--yMin ${mins[@]} \
      #--yMax ${maxs[@]} 
#--sortUsingSamples 1 2 3 \

# plotProfile -m ${outDir}/matrix_${regionData}_${chipData}.mat.gz  \
#               -out ${outDir}/${regionData}_${chipData}_profile.png \
#               --numPlotsPerRow 5 \
#               --regionsLabel "upReg" "downReg" \
#               --startLabel "prom" --endLabel "" \
#               --colors "cyan" "magenta" "brown" "grey" \
#               --yMax ${maxs[@]} \
#               --yMin ${mins[@]} \
#               --samplesLabel ${bwNames[@]} \
# 	      --outFileNameData ${outDir}/${regionData}_${chipData}_profile.tab \
#               --plotTitle "Promoters of COH-1 regulated genes"
