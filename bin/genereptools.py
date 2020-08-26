#!/usr/bin/env python

import sys
import csv
import os.path
import numpy as np
import scipy.stats

# Utils

def file_columns_to_list(filename, c1, c2):
    result = []
    with open(filename, "r") as f:
        for line in f:
            if len(line) > 0 and line[0] != "#":
                line = line.rstrip("\n\r").split("\t")
                result.append([float(line[c1]), float(line[c2])])
    return result

def add_missing(data, maxv):
    # Sort data by descending first element
    data = sorted(data, key=lambda d: d[0], reverse=True)
    result = []
    m = maxv
    for i in range(len(data)):
        tg = data[i][0]
        while m > tg:
            result.append([m, 0])
            m = m - 1
        result.append(data[i])
        m = m - 1
    while m > 0:
        result.append([m, 0])
        m = m - 1
    return result

def conv_and_reverse(data):
    return sorted(data, key=lambda d: d[0], reverse=True)

def safeOpen(files):
    """Open the specified `files' for reading, and return a list of open
streams. Files that do not exist are silently ignored."""
    streams = []
    for f in files:
        if os.path.isfile(f):
            streams.append(open(f, "r"))
    return streams

def compare_mi_histograms(outfile, infile1, infile2, maxv=None):
    """Infile1 = real, infile2 = random."""
    data1 = file_columns_to_list(infile2, 0, 1)
    data2 = file_columns_to_list(infile1, 0, 1)
    # print "{} values read from {}\n{} values read from {}".format(len(data1), infile1, len(data2), infile2)
    # print data1
    # print data2
    tot1 = 0
    tot2 = 0
    maxdiff = [0, 1, 0]
    if maxv:
        data1 = add_missing(data1, maxv)
        data2 = add_missing(data2, maxv)
    else:
        data1 = conv_and_reverse(data1)
        data2 = conv_and_reverse(data2)

    with open(outfile, "w") as out:
        out.write("#Idx\tRandom\tReal\tDiff\tFPR\t% Diff\n")
        for i in range(len(data1)):
            x1 = data1[i][1]
            x2 = data2[i][1]
            tot1 += x1
            tot2 += x2
            diff = tot2-tot1
            # print "{}-{} = {} ({})".format(tot1, tot2, diff, maxdiff)
            if tot2 == 0:
                fpr = 0
            else:
                fpr = 1.0 * tot1 / tot2
            if tot1 == 0:
                pdiff = 0
            else:
                pdiff = 1.0 * diff / tot1
            out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(data1[i][0], tot1, tot2, diff, fpr, pdiff))
            # raw_input()
            if diff > maxdiff[0]:
                maxdiff[0] = diff
                maxdiff[1] = data1[i][0]
                maxdiff[2] = fpr
    return maxdiff

# Commands

class Histograms(object):

    def run(self, args):
        cmd = args[0]
        if cmd == "compare":
            self.runCompare(args[1], args[2], args[3])
        elif cmd == "average":
            self.runAverage(args[1], args[2:])

    def runCompare(self, outfile, realhist, shufhist):
        data1 = file_columns_to_list(shufhist, 0, 1)
        data2 = file_columns_to_list(realhist, 0, 1)
        tot1 = 0
        tot2 = 0
        maxdiff = [0, 1, 0]
        if False: # maxv:
            data1 = add_missing(data1, maxv)
            data2 = add_missing(data2, maxv)
        else:
            data1 = conv_and_reverse(data1)
            data2 = conv_and_reverse(data2)

        with open(outfile, "w") as out:
            out.write("#Idx\tRandom\tReal\tDiff\tFPR\t% Diff\n")
            for i in range(len(data1)):
                x1 = data1[i][1]
                x2 = data2[i][1]
                tot1 += x1
                tot2 += x2
                diff = tot2-tot1
                # print "{}-{} = {} ({})".format(tot1, tot2, diff, maxdiff)
                if tot2 == 0:
                    fpr = 0
                else:
                    fpr = 1.0 * tot1 / tot2
                if tot1 == 0:
                    pdiff = 0
                else:
                    pdiff = 1.0 * diff / tot1
                out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(data1[i][0], tot1, tot2, diff, fpr, pdiff))
                # raw_input()
                if diff > maxdiff[0]:
                    maxdiff[0] = diff
                    maxdiff[1] = data1[i][0]
                    maxdiff[2] = fpr
        return maxdiff

    def runAverage(self, outfile, infiles):
        idx = range(1, 6)
        n = len(infiles)
        streams = safeOpen(infiles)
        if not streams:
            return
        csvs = [ csv.reader(f, delimiter='\t') for f in streams ]
        with open(outfile, "w") as out:
            try:
                for c in csvs:
                    h = c.next()
                out.write("\t".join(h) + "\n")
                while True:
                    try:
                        buf = csvs[0].next()
                    except StopIteration:
                        break
                    for i in idx:
                        buf[i] = float(buf[i])
                    for c in csvs[1:]:
                        buf2 = c.next()
                        for i in idx:
                            buf[i] += float(buf2[i])
                    out.write(buf[0] + "\t" + "\t".join([str(s/n) for s in buf[1:]]) + "\n")
            finally:
                for s in streams:
                    s.close()

class RepeatFilter(object):

    def run(self, args):
        """Remove lines from `source' if they have one value that occurs with frequency higher than `prop'.
Also remove them if their coefficient of variation is less than `mincv'."""
        source = args[0]
        dest = args[1]
        prop = 0.35
        genecols = 2
        mincv = 0.15
        nin = 0
        nbad1 = 0
        nbad2 = 0
        nout = 0
        cvlist = []
        with open(source, "r") as f:
            c = csv.reader(f, delimiter='\t')
            hdr = c.next()
            ncols = len(hdr)-genecols
            maxrepeat = ncols*prop
            with open(dest, "w") as out:
                out.write("\t".join(hdr) + "\n")
                for row in c:
                    nin += 1
                    counts = {}
                    values = []
                    for v in row[genecols:]:
                        x = float(v)
                        if x != 0:
                            values.append(x)
                        if v in counts:
                            counts[v] += 1
                        else:
                            counts[v] = 1
                    nmax = 0
                    for x in counts.itervalues():
                        if x > nmax:
                            nmax = x
                    if nmax < maxrepeat:
                        cv = scipy.stats.variation(values)
                        cvlist.append([row[0], cv])
                        if cv >= mincv:
                            out.write("\t".join(row) + "\n")
                            nout += 1
                        else:
                            nbad2 += 1
                    else:
                        nbad1 += 1
        sys.stdout.write("In={}, Rep={}, CV={}, Out={}\n".format(nin, nbad1, nbad2, nout))

class ZScore(object):
    
    def run(self, infile, outfile):
        with open(outfile, "w") as out:
            with open(infile, "r") as f:
                c = csv.reader(f, delimiter='\t')
                hdr = c.next()
                ncols = len(hdr) - 1
                tail = int(ncols * 0.05)
                out.write("\t".join(hdr) + "\n")
                for line in c:
                    data = [float(v) for v in line[1:]]
                    sdata = sorted(data)[tail:ncols-tail]
                    med = np.median(sdata)
                    std = np.std(sdata)
                    if std == 0:
                        result = [str(v) for v in data]
                    else:
                        result = [ str((v-med)/std) for v in data ]
                    out.write(line[0] + "\t" + "\t".join(result) + "\n")

class SuppFilter(object):
    """Determine optimal support from 'counts' file, generate new consensus."""
    sum = 0                     # Sum of supports
    sumfpr = 0                  # Sum of FPRs
    optsupport = 0
    fprsupport = 0

    def run(self, nrounds, countfiles):
        sumsupp = 0
        sumfpr = 0
        optsupport = 0
        fprsupport = 0
        n = 0
        real = None
        shuffled = []

        for c in countfiles:
            filename = os.path.split(c)[1]
            if filename.startswith("real"):
                real = c
            else:
                shuffled.append(c)
        for sh in shuffled:
            shname = os.path.split(sh)[1].split(".")[0]
            outfile = "real_vs_{}.support.csv".format(shname)
            supp = compare_mi_histograms(outfile, real, sh, maxv=nrounds+1)
            sumsupp += supp[1]
            sumfpr += supp[2]
            n += 1

        optsupport = int(round(1.0 * sumsupp / n))
        fprsupport = 1.0 * sumfpr / n

        with open("optimal-support.csv", "w") as out:
            out.write("OptSupport\t{}\nFPRsupport\t{}\n".format(optsupport, fprsupport))
        
class MIFilter(object):
    
    def run(self, mifiles):
        sum = 0
        sumfpr = 0
        n = 0
        real = None
        shuffled = []

        ## HACK to parse "[a, b, c]" into a proper list of strings.
        ## mifiles = [p.strip() for p in mifiles.strip("[]").split(",")]

        for c in mifiles:
            filename = os.path.split(c)[1]
            if filename.startswith("real"):
                real = c
            else:
                shuffled.append(c)
        for sh in shuffled:
            shname = os.path.split(sh)[1].split(".")[0]
            outfile = "real_vs_{}.mihist.csv".format(shname)
            mi = compare_mi_histograms(outfile, real, sh)
            sum += mi[1]
            sumfpr += mi[2]
            n += 1
        if n > 0:
            optmi = sum / n
            fprmi = sumfpr / n
            with open("optimal-mi.csv", "w") as out:
                out.write("OptMI\t{}\nFPRMI\t{}\n".format(optmi, fprmi))

class SumMIFilter(object):

    def run(self, summifiles):
        sum = 0
        sumfpr = 0
        n = 0
        real = None
        shuffled = []
        
        ## HACK to parse "[a, b, c]" into a proper list of strings.
        ## summifiles = [p.strip() for p in summifiles.strip("[]").split(",")]

        for c in summifiles:
            filename = os.path.split(c)[1]
            if filename.startswith("real"):
                real = c
            else:
                shuffled.append(c)
        for sh in shuffled:
            shname = os.path.split(sh)[1].split(".")[0]
            outfile = "real_vs_{}.summihist.csv".format(shname)
            summi = compare_mi_histograms(outfile, real, sh)
            sum += summi[1]
            sumfpr += summi[2]
            n += 1

        if n > 0:
            optsummi = sum / n
            fprsummi = sumfpr / n
            with open("optimal-summi.csv", "w") as out:
                out.write("OptSumMI\t{}\nFPRSumMI\t{}\n".format(optsummi, fprsummi))

class HistLimits(object):

    def run(self, histfiles):
        histfiles = [p.strip() for p in histfiles.strip("[]").split(",")]
        hmin = np.inf
        hmax = 0
        for hf in histfiles:
            with open(hf, "r") as f:
                c = csv.reader(f, delimiter='\t')
                for line in c:
                    mi = float(line[0])
                    if mi < hmin:
                        hmin = mi
                    if mi > hmax:
                        hmax = mi
        sys.stdout.write("{} {}\n".format(hmin, hmax))

if __name__ == "__main__":
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "histogram":
        Histograms().run(args)
    elif cmd == "repeat":
        RepeatFilter().run(args)
    elif cmd == "zscore":
        ZScore().run(args[0], args[1])
    elif cmd == "suppfilter":
        SuppFilter().run(int(args[0]), args[1:])
    elif cmd == "mifilter":
        MIFilter().run(args)
    elif cmd == "summifilter":
        SumMIFilter().run(args)
    elif cmd == "limits":
        HistLimits().run(args[0])
