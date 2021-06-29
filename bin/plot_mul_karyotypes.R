library(karyoploteR)
library(rtracklayer)
library(GenomeInfoDb)
library(crop)

# Usage Rscript plot_mul_karyotypes.R chm13.bg hg38.bg chm13.sizes hg38.sizes outplot

args = commandArgs(trailingOnly=TRUE)

chm13_sizes = read.table(args[3], header=F, col.names=c('chr', 'len'))
GRCh38_sizes = read.table(args[4], header=F, col.names=c('chr', 'len'))
chroms = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
           'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
           'chr20', 'chr21', 'chr22', 'chrX')
sizes = data.frame(space=chroms, start=rep(0, length(chroms)), end=rep(0, length(chroms)))
c_sizes = data.frame(space=chroms, start=rep(0, length(chroms)), end=rep(0, length(chroms)))
g_sizes = data.frame(space=chroms, start=rep(0, length(chroms)), end=rep(0, length(chroms)))
for (i in 1:length(chroms)){
    sizes$end[i] = max(chm13_sizes$len[chm13_sizes$chr == chroms[i]],
                       GRCh38_sizes$len[GRCh38_sizes$chr == chroms[i]])
    c_sizes$end[i] = chm13_sizes$len[chm13_sizes$chr == chroms[i]]
    g_sizes$end[i] = GRCh38_sizes$len[GRCh38_sizes$chr == chroms[i]]
}
chrom_sizes = toGRanges(sizes, format='other')
c_sizes_gr = toGRanges(c_sizes, format='other')
g_sizes_gr = toGRanges(g_sizes, format='other')
sizes = Seqinfo(seqnames=chroms, seqlengths=sizes$end)

chm13 = read.table(args[1], header=F, col.names=c('chr', 'start', 'end', 'score'))
GRCh38 = read.table(args[2], header=F, col.names=c('chr', 'start', 'end', 'score'))
chm13 = chm13[chm13$score > 0,]
GRCh38 = GRCh38[GRCh38$score > 0,]
chm13$valid = FALSE
GRCh38$valid = FALSE
for (i in 1:length(chroms)){
    chm13$valid[chm13$chr == chroms[i]] = TRUE
    GRCh38$valid[GRCh38$chr == chroms[i]] = TRUE
}
chm13 = chm13[chm13$valid,]
GRCh38 = GRCh38[GRCh38$valid,]
chm13 = subset(chm13, select = -c(valid) )
GRCh38 = subset(GRCh38, select = -c(valid) )

colors = read.table(args[5], header=F, sep="\t", col.names=c('bins', 'cmap'),
                    comment.char="*", colClasses=c('numeric', 'character'))
bins = colors$bins
cmap = colors$cmap
N = length(cmap)
maxval = round(bins[length(bins)])

c_indices = c()
g_indices = c()
for (i in 1:nrow(chm13)){
    c_indices = c(c_indices, findInterval(chm13$score[i], bins) + 1)
}
for (i in 1:nrow(GRCh38)){
    g_indices = c(g_indices, findInterval(GRCh38$score[i], bins) + 1)
}
chm13$color = cmap[c_indices]
GRCh38$color = cmap[g_indices]
write.table(chm13, paste(strsplit(args[1], "\\.bg")[[1]][1], "_color.bg", sep=""),
            sep="\t", row.names=F, col.names=F, quote=F)
write.table(GRCh38, paste(strsplit(args[2], "\\.bg")[[1]][1], "_color.bg", sep=""),
            sep="\t", row.names=F, col.names=F, quote=F)

#for (i in 1:length(chroms)){
#    mask = c()
#    for (j in 1:nrow(chm13)) if (chm13$chr[j] == chroms[i]) mask = c(mask, j)
#    N = length(mask)
#    binsize = mean(chm13$end - chm13$start)
#    chm13$cindex[mask] = chm13$cindex[mask][order(chm13$cindex[mask])]
#    chm13$start[mask] = seq(0,N-1) * binsize + 1
#    for (j in 1:N) chm13$end[mask[j]] = min(c_sizes$end[i], j * binsize + 1)
#    mask = c()
#    for (j in 1:nrow(GRCh38)) if (GRCh38$chr[j] == chroms[i]) mask = c(mask, j)
#    N = length(mask)
#    GRCh38$cindex[mask] = GRCh38$cindex[mask][order(GRCh38$cindex[mask])]
#    GRCh38$start[mask] = seq(0,N-1) * binsize + 1
#    for (j in 1:N) GRCh38$end[mask[j]] = min(c_sizes$end[i], j * binsize + 1)
#}

chm13$strand = "+"
chm13_gr = makeGRangesFromDataFrame(
    chm13,
    seqinfo=sizes,
    seqnames.field='chr',
    start.field='start',
    end.field='end',
    strand.field='strand',
    keep.extra.columns=TRUE)
GRCh38$strand = "-"
GRCh38_gr = makeGRangesFromDataFrame(
    GRCh38,
    seqinfo=sizes,
    seqnames.field='chr',
    start.field='start',
    end.field='end',
    strand.field='strand',
    keep.extra.columns=TRUE)
data = c(chm13_gr, GRCh38_gr)

pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 900
pp$leftmargin = 0.25
pp$ideogramheight = 50
pp$data1height = 400
pp$data2height = 400
pdf(args[6])
kp <- plotKaryotype(genome=chrom_sizes, ideogram.plotter = NULL, plot.type=2,
                    plot.params = pp, chromosomes=chroms)
kpAddCytobandLabels(kp)
kpAddMainTitle(kp, "Mean minimum unique k-mer length per 100Kb", cex=1)
#kpPlotRegions(kp, data=c_sizes_gr, avoid.overlapping = FALSE, col='#bbbbbbff')
#kpPlotRegions(kp, data=g_sizes_gr, avoid.overlapping = FALSE, col='#bbbbbbff',
#              data.panel=2)
kpPlotRegions(kp, data=data[strand(data) == '+'], avoid.overlapping = FALSE,
              col=cmap[c_indices], data.panel=2)
kpPlotRegions(kp, data=data[strand(data) == '-'], avoid.overlapping = FALSE,
              col=cmap[g_indices], data.panel=1)
kpPlotRegions(kp, data=c_sizes_gr, avoid.overlapping = FALSE, col='#ffffff00',
              border='black', data.panel=2, lwd=1)
kpPlotRegions(kp, data=g_sizes_gr, avoid.overlapping = FALSE, col='#ffffff00',
              border='black', data.panel=1, lwd=1)
kpAddLabels(kp, "CHM13", data.panel=2, cex=0.65, col="#000000")
kpAddLabels(kp, "GRCh38", data.panel=1, cex=0.65, col="#000000")
kpAddChromosomeNames(kp, yoffset=-70)
dev.off.crop(args[6])
