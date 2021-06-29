#!/usr/bin/env python3

import sys

import pyBigWig
import numpy

def main():
    bw_fname, out_fname = sys.argv[1:3]
    target = int(sys.argv[3])
    bw = pyBigWig.open(bw_fname)
    chroms = list(bw.chroms().items())
    chroms.sort()
    output = open(out_fname, 'w')
    for chrom, maxlen in chroms:
        print(chrom, file=sys.stderr)
        peaks = []
        n = int(numpy.ceil(maxlen / target))
        binsize = maxlen / n
        for i in range(n):
            start = max(0, int(numpy.round(i * binsize)) - 1)
            end = min(maxlen, int(numpy.round((i + 1) * binsize)) + 1)
            peaks = find_peaks(bw, chrom, start, end, maxlen, target, peaks)
        for s, e in peaks:
            print("{}:{}-{} {}".format(chrom, s, e, e - s), file=output)
    output.close()

def find_peaks(bw, chrom, start, end, maxlen, target, peaks):
    n = min(end - start, 10)
    bins = numpy.round(numpy.linspace(start, end, n + 1)).astype(numpy.int32)
    scores = numpy.array(bw.stats(chrom, start, end, nBins=n, type='max'),
                         dtype=numpy.float32)
    scores[numpy.where(numpy.isnan(scores))] = 0
    where = numpy.where(scores >= target)[0]
    if end - start == n:
        for i in where:
            s = i + start
            if s > 0 and s < maxlen - 1:
                test = numpy.array(bw.values(chrom, s - 1, s + 2))
                test[numpy.where(numpy.isnan(test))] = 0
                if test[1] > test[0] and test[1] > test[2]:
                    peaks.append((s, s + int(test[1]) - 1))
    else:
        for i in where:
            if i == 0:
                if scores[i] >= scores[i + 1]:
                    peaks = find_peaks(bw, chrom, bins[i], bins[i + 1], maxlen,
                                       target, peaks)
            elif i == scores.shape[0] - 1:
                if scores[i] >= scores[i - 1]:
                    peaks = find_peaks(bw, chrom, bins[i], bins[i + 1], maxlen,
                                       target, peaks)
            else:
                if scores[i] >= scores[i - 1] and scores[i] >= scores[i + 1]:
                    peaks = find_peaks(bw, chrom, bins[i], bins[i + 1], maxlen,
                                       target, peaks)
    return peaks

if __name__ == "__main__":
    main()