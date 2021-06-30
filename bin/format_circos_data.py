#!/usr/bin/env python3

import sys

import numpy


def main():
    size_fname, anno_fname, pairs_fname, bg_fname = sys.argv[1:5]
    color_fname, minsize, maxsize, conf_name, out_prefix = sys.argv[5:10]
    minsize, maxsize = int(minsize), int(maxsize)

    # Write colors
    colors = load_colors(color_fname)
    output = open("circos_etc/{}_colors.conf".format(out_prefix), 'w')
    print("<colors>", file=output)
    for i in range(len(colors)):
        print("col_{} = {},{},{}".format(i, int(colors[i][1:3], 16),
                                         int(colors[i][3:5], 16), int(colors[i][5:7], 16)), file=output)
    print("</colors>", file=output)
    output.close()

    # Write ideogram
    chrom_sizes = load_sizes(size_fname)
    chroms = list(chrom_sizes.keys())
    chroms.sort()
    arrays = load_arrays(anno_fname)
    output = open("circos_data/{}_karyotype.txt".format(out_prefix), 'w')
    for chrom in chroms:
        print("chr - hs{} {} {} {} {}".format(
              chrom.lstrip('chr'), chrom.lstrip('chr'), 0,
              chrom_sizes[chrom], chrom), file=output)
    for chrom in chroms:
        where = numpy.where(arrays['chrom'] == chrom.encode('utf8'))[0]
        if where.shape[0] == 0:
            continue
        where2 = where[numpy.where(arrays['end'][where[:-1]] >
                                   arrays['start'][where[1:]])]
        if where2.shape[0] > 0:
            mids = (arrays['end'][where2] + arrays['start'][where2+1]) // 2
            arrays['end'][where2] = mids
            arrays['start'][where2+1] = mids
        if arrays['start'][where[0]] > 0:
            print("band hs{} {} {} 0 {} white".format(
                  chrom.lstrip('chr'), "white", "white",
                  arrays['start'][where[0]]), file=output)
        for i in range(where.shape[0]):
            name = arrays['name'][where[i]].decode('utf8')
            array = arrays['array'][where[i]].decode('utf8')
            start = arrays['start'][where[i]]
            end = arrays['end'][where[i]]
            if i > 0:
                end2 = arrays['end'][where[i - 1]]
                if end2 < start:
                    print("band hs{} {} {} {} {} white".format(
                          chrom.lstrip('chr'), "white", "white", end2,
                          start), file=output)
            print("band hs{} {} {} {} {} {}".format(
                  chrom.lstrip('chr'), name, name, start,
                  end, array), file=output)
        if arrays['end'][where[-1]] < chrom_sizes[chrom]:
            print("band hs{} {} {} {} {} white".format(
                  chrom.lstrip('chr'), "white", "white",
                  arrays['end'][where[-1]], chrom_sizes[chrom]),
                  file=output)
    output.close()

    # Write links
    links = load_data(pairs_fname)
    w = 7500000
    f0 = 0.6
    span = maxsize - minsize
    col_breaks = numpy.round(10**numpy.linspace(
        numpy.log10(minsize), numpy.log10(maxsize),
        len(colors) + 1)[1:]).astype(numpy.int32)
    col_indices = numpy.searchsorted(col_breaks, links['size'])
    output1 = open("circos_data/{}_intra_links.txt".format(out_prefix), 'w')
    output2 = open("circos_data/{}_inter_links.txt".format(out_prefix), 'w')
    for i in range(links.shape[0]):
        c1, s1, e1, m1, c2, s2, e2, m2, size = links[i]
        chrom1 = c1
        chrom2 = c2
        c1, c2 = c1.decode('utf8'), c2.decode('utf8')
        if c1 not in chroms or c2 not in chroms:
            continue
        c1 = c1.replace('chr', 'hs')
        c2 = c2.replace('chr', 'hs')
        f = (size - minsize) / span
        s1 = int(round(m1 - (f * f0 + (1 - f0)) * w))
        e1 = int(round(m1 + (f * f0 + (1 - f0)) * w))
        s2 = int(round(m2 - (f * f0 + (1 - f0)) * w))
        e2 = int(round(m2 + (f * f0 + (1 - f0)) * w))
        if c1 == c2:
            f = output1
        else:
            f = output2
        w1 = numpy.where(arrays['chrom'] == chrom1)[0][0]
        w2 = numpy.where(arrays['chrom'] == chrom2)[0][0]
        if (m1 >= arrays['start'][w1] and m1 < arrays['end'][w1] and
                m2 >= arrays['start'][w2] and m2 < arrays['end'][w2]):
            z = i
        else:
            z = i + links.shape[0]
        col = "col_{}".format(col_indices[i])
        print("{} {} {} {} {} {} color={},z={}".format(
              c1, s1, e1, c2, e2, s2, col, z), file=f)
    output1.close()
    output2.close()

    # Write circos color key
    output = open('circos_data/color_key.txt', 'w')
    for i in range(len(colors)):
        print("{}\t{}".format(col_breaks[i], colors[i]), file=output)
    output.close()

    # Write link peaks
    output = open("circos_data/{}_peaks.txt".format(out_prefix), 'w')
    peaks = []
    w = 250000
    for i in range(links.shape[0]):
        c1, s1, e1, m1, c2, s2, e2, m2, size = links[i]
        s1 = m1 - w // 2
        e1 = m1 + w // 2
        s2 = m2 - w // 2
        e2 = m2 + w // 2
        c1, c2 = c1.decode('utf8'), c2.decode('utf8')
        if c1 not in chroms or c2 not in chroms:
            continue
        size = min(5, numpy.log10(size))
        c1 = c1.replace('chr', 'hs')
        c2 = c2.replace('chr', 'hs')
        col = "col_{}".format(col_indices[i])
        peaks.append((size, c1, s1, e1, col))
        peaks.append((size, c2, s2, e2, col))
    peaks.sort(reverse=True)
    for i in range(len(peaks)):
        size, c, s, e, col = peaks[i]
        size = min(5, size)
        print("{} {} {} {} fill_color={},color={}".format(
            c, s, e, size, col, col), file=output)
    output.close()

    # Write line data
    bg = load_bg(bg_fname, chrom_sizes)
    max_bg = 0
    output = open("circos_data/{}_mu.txt".format(out_prefix), 'w')
    output2 = open("circos_data/{}_missing_mu.txt".format(out_prefix), 'w')
    for chrom in chroms:
        cbg = bg[chrom]
        max_bg = max(max_bg, numpy.amax(cbg['score']))
        c1 = chrom.replace('chr', 'hs')
        print("{} {} {} {}".format(c1, 0, 0, 0), file=output)
        for i in range(cbg.shape[0]):
            s, e, score = cbg[i]
            print("{} {} {} {}".format(c1, s, e, score), file=output)
        print("{} {} {} {}".format(c1, chrom_sizes[chrom],
                                   chrom_sizes[chrom], 0), file=output)
        zeros = numpy.where(cbg['score'] == 0)[0]
        if zeros.shape[0] > 1:
            breaks = numpy.r_[0, numpy.where(zeros[1:] -
                                             zeros[:-1] > 1)[0] + 1,
                              zeros.shape[0]]
            for i in range(breaks.shape[0] - 1):
                s = cbg['start'][zeros[breaks[i]]]
                e = cbg['end'][zeros[breaks[i + 1] - 1]]
                print("{} {} {} 0".format(c1, s, e), file=output2)
        elif zeros.shape[0] == 1:
            s = cbg['start'][zeros[0]]
            e = cbg['end'][zeros[0]]
            print("{} {} {} 0".format(c1, s, e), file=output2)
    max_bg = 4.5
    output.close()
    output2.close()

    # Write conf file
    output = open("circos_etc/{}_circos.conf".format(out_prefix), 'w')
    for line in open(conf_name):
        line = line.rstrip().replace('GENOME', out_prefix).replace("MAXVAL", str(5))
        print(line, file=output)
    output.close()


def load_colors(fname):
    colors = []
    for line in open(fname):
        c = line.rstrip().split()[1]
        colors.append(c)
    return colors


def load_sizes(fname):
    chrom_sizes = {}
    chroms = []
    for line in open(fname):
        line = line.rstrip('\t').split()
        if not line[0].lstrip('chr').isdigit() and line[0] != "chrX":
            continue
        if line[0] == "chrM":
            continue
        chroms.append(line[0])
        chrom_sizes[line[0]] = int(line[1])
    return chrom_sizes


def load_arrays(fname):
    data = []
    for line in open(fname):
        chrom, start, end, fullname = line.rstrip().split()[:4]
        if fullname.count("_arm") > 0:
            continue
        name = fullname.split('(')[0]
        data.append((chrom, int(start), int(end), fullname, name,
                     name.split('_')[0].lower()))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S6'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32),
                                                ('fname', 'S20'),
                                                ('name', 'S20'),
                                                ('array', 'S20')]))
    chroms = numpy.unique(data['chrom'])
    centro = numpy.zeros(chroms.shape[0], dtype=data.dtype)
    for i, chrom in enumerate(chroms):
        where = numpy.where(data['chrom'] == chrom)[0]
        centro[i] = (chrom, numpy.amin(data['start'][where]),
                     numpy.amax(data['end'][where]), 'dred',
                     'dred', 'dred')
        """
        where = numpy.where((data['chrom'] == chrom) &
                            (data['array'] == b'ct'))[0]
        for j in where:
            if (j > 0 and data['chrom'][j - 1] == chrom and
                data['end'][j - 1] > data['start'][j]):
                data['start'][j] = data['end'][j - 1]
            if (j < data.shape[0] - 1 and data['chrom'][j + 1] == chrom and
                data['start'][j + 1] < data['end'][j]):
                data['end'][j] = data['start'][j + 1]
        """
    return centro


def load_data(fname):
    pairs = []
    for line in open(fname):
        line = line.rstrip().split()
        c1, temp = line[0].split(':')
        s1, e1 = temp.split('-')
        s1, e1 = int(s1), int(e1)
        c2, temp = line[2].split(':')
        s2, e2 = temp.split('-')
        s2, e2 = int(s2), int(e2)
        pairs.append((c1, s1, e1, (s1 + e1) // 2, c2, s2, e2, (s2 + e2) // 2,
                      e1 - s1))
    links = numpy.array(pairs, dtype=numpy.dtype([
        ('c1', 'S20'), ('s1', numpy.int32), ('e1', numpy.int32),
        ('m1', numpy.int32), ('c2', 'S20'), ('s2', numpy.int32),
        ('e2', numpy.int32), ('m2', numpy.int32), ('size', numpy.int32)]))
    links = links[numpy.argsort(links['size'])]
    return links


def load_bg(fname, csizes):
    binsize = 200000
    data = {}
    for line in open(fname):
        chrom, start, end, score = line.rstrip().split()
        data.setdefault(chrom, [])
        data[chrom].append((start, end, score))
    for chrom in data.keys():
        if chrom not in csizes:
            continue
        data[chrom] = numpy.array(data[chrom], dtype=numpy.dtype([
            ('start', numpy.int32), ('end', numpy.int32),
            ('score', numpy.float32)]))
        N = int(numpy.ceil((data[chrom]['end'][-1] - 1) / binsize))
        new_data = numpy.zeros(N, dtype=data[chrom].dtype)
        if N > 1:
            bins = numpy.arange(1, N) * binsize
            breaks = numpy.r_[0, numpy.searchsorted(data[chrom]['end'], bins,
                                                    side='right'),
                              data[chrom].shape[0]]
            sizes = breaks[1:] - breaks[:-1]
            new_data['score'] = numpy.bincount(
                numpy.repeat(numpy.arange(N), sizes),
                weights=data[chrom]['score']) / numpy.maximum(1, sizes)
            new_data['start'] = numpy.r_[bins - binsize, bins[-1]]
            new_data['end'] = numpy.r_[bins, data[chrom]['end'][-1]]
        else:
            new_data['score'][0] = numpy.mean(data[chrom]['score'])
            new_data['start'][0] = 0
            new_data['end'][0] = binsize
        new_data['end'][-1] = min(new_data['end'][-1], csizes[chrom])
        data[chrom] = new_data
        data[chrom]['score'] = numpy.log10(
            numpy.maximum(1, data[chrom]['score']))
        data[chrom]['score'] = numpy.minimum(5, data[chrom]['score'])
    return data


if __name__ == "__main__":
    main()
