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
sizes = data.frame(space=c('all'), start=c(0), end=c(0))
c_sizes = data.frame(space=c('all'), start=c(0), end=c(0))
g_sizes = data.frame(space=c('all'), start=c(0), end=c(0))
for (i in 1:length(chroms)){
    c_sizes$end[1] = c_sizes$end[1] + chm13_sizes$len[chm13_sizes$chr == chroms[i]]
    g_sizes$end[1] = g_sizes$end[1] + GRCh38_sizes$len[GRCh38_sizes$chr == chroms[i]]
}
sizes$end[1] = max(c_sizes$end[1] / 10, g_sizes$end[1] / 10)
c_sizes$end[1] = c_sizes$end[1] / 10
g_sizes$end[1] = g_sizes$end[1] / 10
chrom_sizes = toGRanges(sizes, format='other')
c_sizes_gr = toGRanges(c_sizes, format='other')
g_sizes_gr = toGRanges(g_sizes, format='other')
sizes = Seqinfo(seqnames=c("all"), seqlengths=sizes$end)

chm13 = read.table(args[1], header=F,
                   col.names=c('chr', 'size', 'count'))
GRCh38 = read.table(args[2], header=F,
                    col.names=c('chr', 'size', 'count'))
chm13 = chm13[chm13$chr == 'all',]
GRCh38 = GRCh38[GRCh38$chr == 'all',]
chm13$csum = cumsum(chm13$count / 10) 
GRCh38$csum = cumsum(GRCh38$count / 10)

temp = c()
r_a = c(0, 1, 0, 0, -1, 0)
r_b = c(0, 0, 1, 1, 1, 0)
b_a = c(0, 0, -1, 0, 0, 1)
b_b = c(1, 1, 1, 0, 0, 0)
g_a = c(-1, 0, 0, 1, 0, 0)
g_b = c(1, 0, 0, 0, 1, 1)
for (i in 1:6){
  r_a0 = r_a[i]
  r_b0 = r_b[i]
  g_a0 = g_a[i]
  g_b0 = g_b[i]
  b_a0 = b_a[i]
  b_b0 = b_b[i]
  for (j in 1:255){
    r = format(as.hexmode(r_b0 * 255 + r_a0 * j), width=2, upper.case=T)
    g = format(as.hexmode(g_b0 * 255 + g_a0 * j), width=2, upper.case=T)
    b = format(as.hexmode(b_b0 * 255 + b_a0 * j), width=2, upper.case=T)
    temp = c(temp, paste("#", r, g, b, "FF", sep=""))
  }
}
temp = c(temp[200:length(temp)], temp[1:199])
temp = temp[1:round(length(temp) * 0.9)]
cmap = c()
pos = 1
for (i in 1:3){
  for (j in 1:(2^(4 + i))){
    cmap = c(cmap, temp[pos])
    pos = pos + 2^(4 - i)
  }
}
cmap = c(cmap, temp[pos:length(temp)])

minv = 17
maxv = max(c(chm13$size, GRCh38$size))
N = length(cmap)
n = N - 240 - 1
m = minv + 239
step = (log10(maxv / m) - log10((m + 1) / m)) / n
bins = c(minv:(minv + 239), m * 10^((step * 0:n) + log10((m + 1) / m)))
key = data.frame(bin=bins, color=cmap)
fname = paste(strsplit(args[5], "\\.pdf")[[1]][1], "_key.txt", sep="")
write.table(key, fname, sep="\t", row.names=F, col.names=F, quote=F)

c_data = data.frame(chr=rep('all', N), start=rep(0, N), end=rep(0, N),
                    strand=rep('+', N), index=(1:N))
g_data = data.frame(chr=rep('all', N), start=rep(0, N), end=rep(0, N),
                    strand=rep('-', N), index=(1:N))
for (i in 1:(N - 1)){
    c_data$end[i] = sum(chm13$count[chm13$size <= bins[i]] / 10)
    c_data$start[i+1] = c_data$end[i]
    g_data$end[i] = sum(GRCh38$count[GRCh38$size <= bins[i]] / 10)
    g_data$start[i+1] = g_data$end[i]
}
i = 1
while (is.nan(c_data$end[i])){
    c_data$start[i] = 0
    c_data$end[i] = 0
    i = i + 1
}
c_data$start[i] = 0
i = 1
while (is.nan(g_data$end[i])){
    g_data$start[i] = 0
    g_data$end[i] = 0
    i = i + 1
}
g_data$start[i] = 0
c_data$end[nrow(c_data)] = sum(chm13$count / 10)
g_data$end[nrow(g_data)] = sum(GRCh38$count / 10)
c_data = c_data[c_data$end > c_data$start,]
g_data = g_data[g_data$end > g_data$start,]

kmers = c(15, 20, 25, 35, 50, 100, 300, 1000, 10000, 97700)
c_mids = c()
c_labels = c("15", "20", "25", "35", "50", "100", "300", "1000", "10000", "9.77e4")
for (i in 1:length(kmers)){
    k = kmers[i]
    if (k <= nrow(chm13)){
        c_mids = c(c_mids, sum(chm13$count[1:(k-1)] / 10) + chm13$count[k] / 20)
    }
}
kmers = c(15, 20, 25, 35, 50, 100, 300, 1000, 10000, 113000)
g_mids = c()
g_labels = c()
g_labels = c("15", "20", "25", "35", "50", "100", "300", "1000", "10000", "1.13e5")
for (i in 1:length(kmers)){
    k = kmers[i]
    if (k <= nrow(GRCh38)){
        g_mids = c(g_mids, sum(GRCh38$count[1:(k-1)] / 10) + GRCh38$count[k] / 20)
    }
}

data = makeGRangesFromDataFrame(
    rbind(c_data, g_data),
    seqinfo=sizes,
    seqnames.field='chr',
    start.field='start',
    end.field='end',
    strand.field='strand',
    keep.extra.columns=TRUE)
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 01
pp$data1inmargin <- 01
pp$data2outmargin <- 01
pp$data2inmargin <- 01
pp$topmargin <- 45
pp$bottommargin <- 150
pp$leftmargin = 0.24
pp$ideogramheight = 1
pp$data1height = 15
pp$data2height = 15
pdf(args[5])
kp <- plotKaryotype(genome=chrom_sizes, plot.type=2, ideogram.plotter=NULL,
                    plot.params=pp, main="Minimum unique k-mer length (bp)",
                    cex=1.5, labels.plotter=NULL)
kpPlotRegions(kp, data=data[strand(data) == "+"], avoid.overlapping = FALSE,
              col=cmap[c_data$index])
kpPlotRegions(kp, data=data[strand(data) == "-"], data.panel=2,
              avoid.overlapping = FALSE, col=cmap[g_data$index])
kpPlotRegions(kp, data=c_sizes_gr, avoid.overlapping = FALSE, col='#ffffff00',
              border='black', lwd=1)
kpPlotRegions(kp, data=g_sizes_gr, avoid.overlapping = FALSE, col='#ffffff00',
              border='black', data.panel=2, lwd=1)
kpAddLabels(kp, "CHM13", cex=1.5, col="black")
kpAddLabels(kp, "GRCh38", data.panel=2, cex=1.5, col="black")
kpPlotMarkers(kp, chr=rep('all', length(c_mids)), x=c_mids, labels=c_labels,
              text.orientation='vertical', r0=1, r1=1.5, cex=1.25,
              marker.parts=c(0.75,0.75,0.4), label.margin=0.2, pos=2, srt=270,
              label.dist=0.002, ignore.chromosome.ends=TRUE)
kpPlotMarkers(kp, chr=rep('all', length(g_mids)), x=g_mids, labels=g_labels,
              text.orientation='vertical', data.panel=2, r0=1.0, r1=1.5,
              cex=1.25, marker.parts=c(0.75,0.75,0.4), label.margin=0.2,
              pos=4, srt=270, label.dist=0.002, ignore.chromosome.ends=TRUE)
dev.off.crop(args[5])
