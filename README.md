# FischleLab_Eswara_ChIP

Scripts for analysis of ChIP seq. Additional data for:
>HP1 dimerization but not LLPS drives the condensation and segregation of H3K9me-marked chromatin. Karthik Eswara, Jennifer Semple, Francesca Rivas-Cuestas, Shakila Ali, Sonia El Mouridi, Yogesh Ostwal, Christian Frøkjær-Jensen, Peter Meister, Wolfgang Fischle
doi: https://doi.org/10.1101/2025.11.13.687974

1. ***makeSampleSheet.R*** - script to create input sample sheet for nf-core pipeline

2. ***nf_chipseq_submit.sh*** - used to run [nf-core/chipseq](https://nf-co.re/chipseq/2.1.0/) pipeline version 2.1.0

3. ***gatherQCplots.py*** - to collect individual pdfs made by nf-core into a single file

4. ***averageBigwigs.sh*** - script to average repeats and create ip-input and log2(ip/input) tracks, using deeptools.

5. ***prepareForDeeptools.R*** - script to create a file with a list of average log2(IP/input) with separate column for sample name to use as input for deeptools. Also it takes the table of consensus peaks output by nf-core/chipseq pipeline and creates a bedfile of peaks present at least 3 replicates among all H3K9me2 ChIP seq samples.

6. ***H3K9me2_overEM88peaks.sh*** - deeptools showing H3K9me2 signal from different samples over set of consensus peaks.
