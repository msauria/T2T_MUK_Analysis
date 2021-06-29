library(karyoploteR)
library(rtracklayer)
library(GenomeInfoDb)
library(crop)

#### Plot Key ####

args = commandArgs(trailingOnly=TRUE)

colors = read.table(args[1], header=F, sep="\t", col.names=c('bins', 'cmap'),
                    comment.char="*", colClasses=c('numeric', 'character'))
bins = colors$bins
cmap = colors$cmap
N = length(cmap)

sizes = data.frame(space=c('key'), start=c(0), end=c(N * 10000))
chrom_sizes = toGRanges(sizes, format='other')
sizes = Seqinfo(seqnames=c('key'), seqlengths=c(N * 10000))
pos = seq(1.5, 6, 0.5)
n = length(pos)
cyto = data.frame(chr=rep('key', N), start=((0:(N-1)) * 10000),
                  end=((1:N) * 10000),
                  name=as.character(rep('', N)), gieStain=rep('geng', N))
cyto_gr = makeGRangesFromDataFrame(
    cyto,
    seqinfo=sizes,
    seqnames.field='chr',
    start.field='start',
    end.field='end',
    keep.extra.columns=TRUE)

kmers = c(20, 50, 100, 150, 300, 1E3, 1E4, 1E5)
k_labels = c("20", "50", "100", "150", "300", "1e3", "1e4", "1e5")
k_pos = c()
for (i in 1:length(kmers)){
    k_pos = c(k_pos, ((N - findInterval(kmers[i], bins) + 0.5)) * 10000)
}

pp <- getDefaultPlotParams(plot.type=4)
pp$ideogramheight=0.5
pp$leftmargin=0.05
pp$rightmargin=0.05
pp$data1height=1.25
pp$topmargin=2
pdf(args[2])
kp <- plotKaryotype(genome=chrom_sizes, ideogram.plotter=NULL, plot.type=4,
                    plot.params = pp, labels.plotter=NULL)
kpAddMainTitle(kp, "Mean minimum unique k-mer length (bp)", cex=1.5)
kpPlotRegions(kp, data=cyto_gr, avoid.overlapping = FALSE, col=cmap[N:1],
              data.panel=2)
kpPlotRegions(kp, data=chrom_sizes, avoid.overlapping = FALSE,
              col="#ffffff00", border="#000000")
kpPlotMarkers(kp, chr=rep('key', length(k_pos)), x=k_pos, y=1,
              labels=k_labels, r0=0, r1=-0.5, pos=2,
              marker.parts=c(0.1,0.1,0.1), label.margin=-0.2, cex=1.5,
              ignore.chromosome.ends=TRUE)
dev.off.crop(args[2])
