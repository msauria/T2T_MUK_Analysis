#!/usr/bin/env python3

import sys
import struct

def main():
    fasta_fname, sa_fname, target_fname, out_fname, minsize = sys.argv[1:6]
    minsize = int(minsize)
    fa, chr_indices, chroms = load_fasta(fasta_fname)
    sa = open(sa_fname, 'rb')
    pairs = set()
    N = chr_indices[-1]
    output = open(out_fname, 'w')
    for line in open(target_fname):
        target, span = line.rstrip().split()
        chrom, temp = target.split(":")
        start, end = temp.split("-")
        start, end = int(start), int(end)
        span = end - start
        index = start + chr_indices[chroms.index(chrom)]
        seq = fa[index:(index + span)]
        p, s, e = 0, 0, 2*N
        iteration = 0
        while e - s > 1:
            p, s, e = find_sa_index(sa, p, s, e, seq, fa)
            iteration += 1
            if iteration >= 40:
                break
        sa.seek((s - p - 1)*8, 1)
        sa_index0 = struct.unpack("Q", sa.read(8))[0]
        sa_index1 = struct.unpack("Q", sa.read(8))[0]
        sa_index2 = struct.unpack("Q", sa.read(8))[0]
        chrom0, index0, seq0 = get_seq(fa, sa_index0, span, chr_indices, chroms)
        chrom1, index1, seq1 = get_seq(fa, sa_index1, span, chr_indices, chroms)
        chrom2, index2, seq2 = get_seq(fa, sa_index2, span, chr_indices, chroms)
        if sa_index0 < N:
            strand0 = "for"
        else:
            strand0 = "rev"
        if sa_index1 < N:
            strand1 = "for"
        else:
            strand1 = "rev"
        if sa_index2 < N:
            strand2 = "for"
        else:
            strand2 = "rev"
        lcp0 = 0
        while lcp0 < len(seq0) and seq0[lcp0] == seq1[lcp0]:
            lcp0 += 1
        lcp1 = 0
        while lcp1 < len(seq1) and seq1[lcp1] == seq2[lcp1]:
            lcp1 += 1
        if lcp0 >= minsize:
            if (chrom0, index0) > (chrom1, index1):
                chrom0, chrom1 = chrom1, chrom0
                index0, index1 = index1, index0
                strand0, strand1 = strand1, strand0
            key = (chrom0, index0, chrom1, index1)
            if key not in pairs:
                print("{}:{}-{} {}\t{}:{}-{} {}".format(
                    chrom0, index0, index0+lcp0, strand0,
                    chrom1, index1, index1+lcp0, strand1), file=output)
                pairs.add(key)
        if lcp1 >= minsize:
            if (chrom1, index1) > (chrom2, index2):
                chrom1, chrom2 = chrom2, chrom1
                index1, index2 = index2, index1
                strand1, strand2 = strand2, strand1
            key = (chrom1, index1, chrom2, index2)
            if key not in pairs:
                print("{}:{}-{} {}\t{}:{}-{} {}".format(
                    chrom1, index1, index1+lcp1, strand1,
                    chrom2, index2, index2+lcp1, strand2), file=output)
                pairs.add(key)
        sa.seek(0, 0)
    output.close()

def load_fasta(fname):
    chr_indices = []
    chroms = []
    fa = []
    pos = 0
    for line in open(fname):
        if line.startswith('>'):
            chroms.append(line[1:].rstrip())
            chr_indices.append(pos)
        else:
            fa.append(line.rstrip().upper())
            pos += len(fa[-1])
    chr_indices.append(pos)
    return "".join(fa), chr_indices, chroms

def find_sa_index(sa, pos, start, end, seq, fa):
    mid = (end + start) // 2
    sa.seek((mid - pos)*8, 1)
    index = struct.unpack("Q", sa.read(8))[0]
    pos = mid + 1
    if index >= len(fa):
        tseq = revcomp(fa, index, len(seq))
    else:
        tseq = fa[index:min(index + len(seq), len(fa))]
    if seq > tseq:
        return pos, mid + 1, end
    elif seq < tseq:
        return pos, start, mid
    else:
        return pos, mid, mid+1

def revcomp(fa, index, span):
    if index > len(fa):
        index2 = len(fa) * 2 - index - span
    else:
        index2 = index
    seq = []
    for i in range(max(0, index2), index2 + span)[::-1]:
        base = fa[i]
        seq.append(convert_base(base))
    return "".join(seq)

def convert_base(base):
    if base == "A":
        return "T"
    elif base == "T":
        return "A"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    else:
        return "N"

def get_seq(fa, index, span, chr_indices, chroms):
    if index < len(fa):
        index2 = index
        seq = fa[index:(index + span)]
    else:
        index2 = len(fa)*2 - index - span
        seq = revcomp(fa, index, span)
    chrint = 0
    while index2 > chr_indices[chrint+1]:
        chrint += 1
    return (chroms[chrint], index2 - chr_indices[chrint], seq)



if __name__ == "__main__":
    main()