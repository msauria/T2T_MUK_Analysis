#!/usr/bin/env python3

import sys

import pyBigWig
import numpy

def main():
    mu_fname, out_fname, binsize = sys.argv[1:4]
    binsize = int(binsize)
    mu = pyBigWig.open(mu_fname)
    chroms = list(mu.chroms().items())
    chroms.sort()
    output = open(out_fname, 'w')
    scores = {}
    for chrom, length in chroms:
        n = int(numpy.round((length - 1) / binsize)) + 1
        score = numpy.array(mu.stats(chrom, 0, length - 1, type='mean',
                            exact=True, nBins=n), numpy.float32)
        score[numpy.where(numpy.isnan(score))] = 0
        for i, s in enumerate(score):
            if s > 0:
                print("{}\t{}\t{}\t{}".format(chrom, binsize*i,
                                              binsize*(i+1), s), file=output)
    output.close()

if __name__ == "__main__":
    main()
