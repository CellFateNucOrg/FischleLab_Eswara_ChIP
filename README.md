# FischleLab_Eswara_ChIP

Scripts for analysis of ChIP seq.

1. ***makeSampleSheet.R*** script to create input sample sheet for nf-core pipeline
2. ***nf_chipseq_submit.sh*** used to run nf-core pipeline
3. ***gatherQCplots.py*** to collect individual pdfs made by nf-core into a single file
4. ***averageBigwigs.sh*** script to average repeats and create ip-input and log2(ip/input) tracks, using deeptools
5. ***subsetChIPpeaks.R*** script to create subsets of nf-core's consensus peaks to determine best peaks to use.
