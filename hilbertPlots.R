#library(HilbertVis)
#library(HilbertVisGUI)
library(GenomeInfoDb)
library(BSgenome.Celegans.UCSC.ce11)
library(HilbertCurve)
library(ComplexHeatmap)
library(rtracklayer)
library(circlize)

n2<-import("/Volumes/meister.data/FischleLab_KarthikEswara/ChIPseq_H3K9me2_20260405/bigwigLog2fc/N2_IPlog2fcInput_avr.bigwig")
n2<-convertToCe11(n2)

EM88<-import("/Volumes/meister.data/FischleLab_KarthikEswara/ChIPseq_H3K9me2_20260405/bigwigLog2fc/EM88_IPlog2fcInput_avr.bigwig")
EM88<-convertToCe11(EM88)

EM90<-import("/Volumes/meister.data/FischleLab_KarthikEswara/ChIPseq_H3K9me2_20260405/bigwigLog2fc/EM90_IPlog2fcInput_avr.bigwig")
EM90<-convertToCe11(EM90)

EM91<-import("/Volumes/meister.data/FischleLab_KarthikEswara/ChIPseq_H3K9me2_20260405/bigwigLog2fc/EM91_IPlog2fcInput_avr.bigwig")
EM91<-convertToCe11(EM91)

EM92<-import("/Volumes/meister.data/FischleLab_KarthikEswara/ChIPseq_H3K9me2_20260405/bigwigLog2fc/EM92_IPlog2fcInput_avr.bigwig")
EM92<-convertToCe11(EM92)

convertToCe11<-function(gr){
  seqlevelsStyle(gr)<-"UCSC"
  seqlevels(gr)<-seqlevels(Celegans)
  sort(gr)
  return(gr)
}



drawHilbertCurve<-function(gr,chr,title=NULL,legendTitle=NULL, col_fun=NULL,lgd=NULL){
  if(is.null(col_fun)){
    col_fun = colorRamp2(quantile(gr$score, c(0.1, 0.5, 0.9)), c("blue", "#FFFFCC", "red"))
    lgd = Legend(col_fun = col_fun, title = legendTitle,
               at = quantile(gr$score, c(0, 0.5, 1)),
               labels = round(quantile(gr$score, c(0,  0.5,  1)),2))
  }
  hc=GenomicHilbertCurve(chr=chr,level=9,species="ce11",mode="pixel",legend=lgd,
                         title=paste(title,chr))
  hc_layer(hc, gr, col = col_fun(gr$score))
}



col_fun = colorRamp2(quantile(n2$score, c(0.1, 0.5, 0.9)), c("blue", "#FFFFCC", "red"))
lgd = Legend(col_fun = col_fun, title = "H3K9me2",
             at = quantile(n2$score, c(0, 0.5, 1)),
             labels = round(quantile(n2$score, c(0,  0.5,  1)),2))



pdf("hilbertCurves_chrII_sameScale.pdf",width=8,height=11)
drawHilbertCurve(n2,chr="chrII",title="N2",legendTitle="H3K9me2",col_fun,lgd)
drawHilbertCurve(EM88,chr="chrII",title="EM88",legendTitle="H3K9me2",col_fun,lgd)
drawHilbertCurve(EM90,chr="chrII",title="EM90",legendTitle="H3K9me2",col_fun,lgd)
drawHilbertCurve(EM91,chr="chrII",title="EM91",legendTitle="H3K9me2",col_fun,lgd)
drawHilbertCurve(EM92,chr="chrII",title="EM92",legendTitle="H3K9me2",col_fun,lgd)
dev.off()

# hc=GenomicHilbertCurve(chr="chrI",level=6,species="ce11",legend=lgd)
# hc_points(hc,n2,gp = gpar(fill = col_fun(n2$score), col = col_fun(n2$score)))
#
# hc=GenomicHilbertCurve(chr="chrI",level=9,species="ce11",mode="pixel",legend=lgd,
#                        title="N2 chrI")
# hc_layer(hc,n2, col = col_fun(n2$score))



# cov<-coverage(n2,width=1000,weight="score")
#
# plotLongVector(cov)
# chrI<-n2[seqnames(n2)=="I"]$score
# plotLongVector(chrI)
# plotHilbertCurve(chrI)
#
# hilbertDisplay(cov$"I")
#
#
# library( grid )
# pushViewport( viewport( layout=grid.layout( 2, 2 ) ) )
# for( i in 1:4 ) {
#   pushViewport( viewport(
#     layout.pos.row=1+(i-1)%/%2, layout.pos.col=1+(i-1)%%2 ) )
#   plotHilbertCurve( i, new.page=FALSE )
#   popViewport( )
# }
