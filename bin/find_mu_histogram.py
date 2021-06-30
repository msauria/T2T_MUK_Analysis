#!/usr/bin/env python3

import sys

import pyBigWig
import numpy


def main():
    mu_fname, out_fname = sys.argv[1:3]
    mu = pyBigWig.open(mu_fname)
    chroms = list(mu.chroms().items())
    chroms.sort()
    allcounts = numpy.zeros(mu.header()['maxVal'] + 1, dtype=numpy.int64)
    output = open(out_fname, 'w')
    for chrom, length in chroms:
        # if not chrom.lstrip('chr').isdigit() and chrom != "chrX":
        #    continue
        counts = numpy.zeros(mu.header()['maxVal'] + 1, dtype=numpy.int64)
        print(chrom)
        scores = numpy.array(mu.values(chrom, 0, length), numpy.float64)
        scores = scores[numpy.where(numpy.logical_not(numpy.isnan(scores)))]
        counts += numpy.bincount(scores.astype(numpy.int64),
                                 minlength=counts.shape[0])
        allcounts += counts
        for i in range(counts.shape[0]):
            print("{}\t{}\t{}".format(chrom, i + 1, counts[i]), file=output)
    for i in range(allcounts.shape[0]):
        print("all\t{}\t{}".format(i + 1, allcounts[i]), file=output)
    output.close()


if __name__ == "__main__":
    main()
