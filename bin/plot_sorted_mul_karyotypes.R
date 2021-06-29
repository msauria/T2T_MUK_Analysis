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

colors = read.table(args[5], header=F, sep="\t", col.names=c('bins', 'cmap'),
                    comment.char="*", colClasses=c('numeric', 'character'))
bins = colors$bins
cmap = colors$cmap
N = length(cmap)
maxval = round(bins[length(bins)])
chm13 = read.table(args[1], header=F, col.names=c('chr', 'size', 'count'))
GRCh38 = read.table(args[2], header=F, col.names=c('chr', 'size', 'count'))
chm13 = chm13[chm13$chr != 'all' & chm13$count > 0,]
GRCh38 = GRCh38[GRCh38$chr != 'all' & GRCh38$count > 0,]

chm13$bindex = 0
for (i in 1:length(chroms)){
    print(chroms[i])
    cmask = chm13$chr == chroms[i]
    ctemp = (1:nrow(chm13))[cmask]
    cstart = ctemp[1]
    cend = ctemp[length(ctemp)]
    cindices = c(cstart - 1)
    for (j in 1:N){
        cindices = c(cindices, cstart + findInterval(bins[j], chm13$size[cstart:cend]) - 1)
        if (cindices[j + 1] > cindices[j]){
            chm13$bindex[(cindices[j] + 1):cindices[j + 1]] = j
        }
    }
}

GRCh38$bindex = 0
for (i in 1:length(chroms)){
    print(chroms[i])
    gmask = GRCh38$chr == chroms[i]
    gtemp = (1:nrow(GRCh38))[gmask]
    gstart = gtemp[1]
    gend = gtemp[length(gtemp)]
    gindices = c(gstart - 1)
    for (j in 1:N){
        gindices = c(gindices, findInterval(bins[j], GRCh38$size[gstart:gend]) + gstart - 1)
        if (gindices[j + 1] > gindices[j]){
            GRCh38$bindex[(gindices[j] + 1):gindices[j + 1]] = j
        }
    }
}

chm13_b = data.frame(chr=rep(chroms, each=N),
                     start=rep(0, N*length(chroms)),
                     end=rep(0, N*length(chroms)),
                     cindex=rep(1:N, length(chroms)))
GRCh38_b = data.frame(chr=rep(chroms, each=N),
                      start=rep(0, N*length(chroms)),
                      end=rep(0, N*length(chroms)),
                      cindex=rep(1:N, length(chroms)))
for (i in 1:length(chroms)){
    print(chroms[i])
    cmask = chm13$chr == chroms[i]
    ctemp = (1:nrow(chm13))[cmask]
    cstart = ctemp[1]
    cend = ctemp[length(ctemp)]
    cmask = chm13$bindex[cstart:cend] == 1
    ctemp = (cstart:cend)[cmask]
    total = 0
    if (length(ctemp) > 0){
        cstart1 = ctemp[1]
        cend1 = ctemp[length(ctemp)]
        total = sum(chm13$count[cstart1:cend1])
        chm13_b$end[(i - 1) * N + 1] = sum(chm13$count[cstart1:cend1])
    } else {
        cstart1 = cstart
        cend1 = cstart
    }
    gmask = GRCh38$chr == chroms[i]
    gtemp = (1:nrow(GRCh38))[gmask]
    gstart = gtemp[1]
    gend = gtemp[length(gtemp)]
    gmask = GRCh38$bindex[gstart:gend] == 1
    gtemp = (gstart:gend)[gmask]
    if (length(gtemp) > 0){
        gstart1 = gtemp[1]
        gend1 = gtemp[length(gtemp)]
        GRCh38_b$end[(i - 1) * N + 1] = sum(GRCh38$count[gstart1:gend1])
    }
    for (j in 2:N){
        index1 = (i - 1) * N + j
        index0 = index1 - 1
        cmask = chm13$bindex[cstart:cend] == j
        ctemp = (cstart:cend)[cmask]
        gmask = GRCh38$bindex[gstart:gend] == j
        gtemp = (gstart:gend)[gmask]
        chm13_b$start[index1] = chm13_b$end[index0]
        if (length(ctemp) > 0){
            cstart1 = ctemp[1]
            cend1 = ctemp[length(ctemp)]
            total = total + sum(chm13$count[cstart1:cend1])
            chm13_b$end[index1] = (sum(chm13$count[cstart1:cend1])
                + chm13_b$end[index0])
        } else chm13_b$end[index1] = chm13_b$end[index0]
        GRCh38_b$start[index1] = GRCh38_b$end[index0]
        if (length(gtemp) > 0){
            gstart1 = gtemp[1]
            gend1 = gtemp[length(gtemp)]
            GRCh38_b$end[index1] = (sum(GRCh38$count[gstart1:gend1])
                + GRCh38_b$end[index0])
        } else GRCh38_b$end[index1] = GRCh38_b$end[index0]
    }
}
chm13_b = chm13_b[chm13_b$end > chm13_b$start,]
GRCh38_b = GRCh38_b[GRCh38_b$end > GRCh38_b$start,]

chm13_b$strand = "+"
chm13_gr = makeGRangesFromDataFrame(
    chm13_b,
    seqinfo=sizes,
    seqnames.field='chr',
    start.field='start',
    end.field='end',
    strand.field='strand',
    keep.extra.columns=TRUE)
GRCh38_b$strand = "-"
GRCh38_gr = makeGRangesFromDataFrame(
    GRCh38_b,
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
kpAddMainTitle(kp, "Minimum unique k-mer length (bp)", cex=1)
kpPlotRegions(kp, data=data[strand(data) == '+'], avoid.overlapping = FALSE,
              col=cmap[chm13_b$cindex], data.panel=2)
kpPlotRegions(kp, data=data[strand(data) == '-'], avoid.overlapping = FALSE,
              col=cmap[GRCh38_b$cindex], data.panel=1)
kpPlotRegions(kp, data=c_sizes_gr, avoid.overlapping = FALSE, col='#ffffff00',
              border='black', data.panel=2, lwd=1)
kpPlotRegions(kp, data=g_sizes_gr, avoid.overlapping = FALSE, col='#ffffff00',
              border='black', data.panel=1, lwd=1)
kpAddLabels(kp, "CHM13", data.panel=2, cex=0.65, col="#000000")
kpAddLabels(kp, "GRCh38", data.panel=1, cex=0.65, col="#000000")
kpAddChromosomeNames(kp, yoffset=-70)
dev.off.crop(args[6])
