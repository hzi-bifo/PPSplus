#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.


    Summary visualization of the simulated datasets.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from algbioi.com import fasta as fas
from algbioi.com import csv
from algbioi.com import taxonomy_ncbi


def getData():
    """
        @return: data[(method, rank, correction, simset)] = ((cTp/cT)*100, (cW/cT)*100, (cU/cT)*100)
    """
    srcDir = '/Users/ivan/Documents/work/Doc/PPSplus/evaluation_simset/summary/src'
    ranks = ['species', 'genus', 'family']
    LOGNORM_SIZE = 68751.268
    UNIFORM_SIZE = 142570.869

    data = {}
    predFile = None
    for rank in ranks:

        for f in os.listdir(srcDir):
            cTp = 0.0  # cumulative: true predicted
            cW = 0.0   # cumulative: wrong
            cU = 0.0   # cumulative: unassigned
            cT = 0.0   # cumulative: true
            method = f.split('_')[0]
            simset = f.split('_')[3].split('.')[0]
            correction = f.split('.')[-2].split('_')[-1]
            if correction != 'correction':
                correction = 'nc'
            #print method, simset, correction
            if simset == 'lognorm':  # lognorm size
                t = LOGNORM_SIZE
            elif simset == 'uniform':  # ?? uniform size
                t = UNIFORM_SIZE
            else:
                raise TypeError('wrong set type %s' % simset)

            coef = 1.0 # the weight of the current entry
            for line in open(os.path.join(srcDir, f)):
                if line.startswith('@'):
                    coef = float(line.split(',')[1])
                elif line.startswith('$'):
                    predFile = line.strip().split(':')[1]
                elif line.startswith(rank):
                    pr = float(line.split(',')[-6]) / 100.0
                    re = float(line.split(',')[-5]) / 100.0
                    tp = re * t
                    #
                    # p = predBp / 1000.
                    if pr > 1.1 and re > 1.1:
                        p = tp / pr  # TODO: !!!
                    else:
                        print predFile,
                        predBp = assignedAtRank(rank, predFile, simset)
                        p = predBp / 1000.
                        print ' ..ok'
                        if re < 0.05 or pr < 0.05:
                            print rank, p, re, pr


                    w = p - tp
                    u = t - p
                    #print tp, w, u

                    cTp += coef * tp
                    cW += coef * w
                    cU += coef * u
                    cT += coef * t

            if cT > 0.0:
                #print method, rank, correction, simset, (cTp/cT), (cW/cT), (cU/cT)
                data[(method, rank, correction, simset)] = ((cTp/cT)*100, (cW/cT)*100, (cU/cT)*100)
            assert int(cTp + cW + cU) == int(cT)

    return data


def createChart():
    """
        Creates a summary figure for simulated datasets divided into 4 subplots.
    """
    # outFile = '/Users/ivan/Documents/work/Doc/PPSplus/revision/figures_simset/Fig_1_simset_summary_Kraken.png'
    outFile = '/Users/ivan/Documents/work/Doc/PPSplus/peerj/revision1/fig/Fig_1_simset_summary_peerj.png'
    dpi = 300

    N = 3  # N ranks
    ind = np.arange(N)    # the x locations for the groups
    width = 0.16  # 0.2
    width2 = 0.18  # 0.22
    data = getData()  # (method, rank, correction, simset) # ((cTp/cT), (cW/cT), (cU/cT)) true, wrong, unassigned

    # create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, subplot_kw={})  # subplot_kw ... arg passed to add_subplot

    # figure dimensions
    fig.set_figwidth(8)
    fig.set_figheight(8)

    # padding/margins
    plt.subplots_adjust(left=0.1, bottom=0.22, right=0.9, top=0.9, wspace=0.3, hspace=0.75)  # left=0.1, bottom=0.22, right=0.9, top=0.9, wspace=0.3, hspace=0.75

    # create 4 subplots
    for simset, correction, ax, title, letter in [('uniform', 'nc', ax1, 'Uniform dataset', 'A)'),
                                                  ('uniform', 'correction', ax2, 'Uniform dataset, correction', 'B)'),
                                                  ('lognorm', 'nc', ax3, 'Log-norm dataset', 'C)'),
                                                  ('lognorm', 'correction', ax4, 'Log-norm dataset, correction', 'D)')]:

        # add zero data for pps
        if ('pps', 'species', correction, simset) not in data:
            data[('pps', 'species', correction, simset)] = (0.0, 0.0, 100.0)

        # PPS+ (species, genus, family)
        u_ppsp = (100.0, 100.0, 100.0)
        w_ppsp = (data[('ppsp', 'species', correction, simset)][1] + data[('ppsp', 'species', correction, simset)][0],
                  data[('ppsp', 'genus', correction, simset)][1] + data[('ppsp', 'genus', correction, simset)][0],
                  data[('ppsp', 'family', correction, simset)][1] + data[('ppsp', 'family', correction, simset)][0])  # false
        tp_ppsp = (data[('ppsp', 'species', correction, simset)][0], data[('ppsp', 'genus', correction, simset)][0],
                   data[('ppsp', 'family', correction, simset)][0])  # true

        # PPS (species, genus, family)
        u_pps = (100.0, 100.0, 100.0)
        w_pps = (data[('pps', 'species', correction, simset)][1] + data[('pps', 'species', correction, simset)][0],
                  data[('pps', 'genus', correction, simset)][1] + data[('pps', 'genus', correction, simset)][0],
                  data[('pps', 'family', correction, simset)][1] + data[('pps', 'family', correction, simset)][0])  # false
        tp_pps = (data[('pps', 'species', correction, simset)][0], data[('pps', 'genus', correction, simset)][0],
                  data[('pps', 'family', correction, simset)][0])  # true

        # Megan (species, genus, family)
        u_megan = (100.0, 100.0, 100.0)
        w_megan = (data[('megan', 'species', correction, simset)][1] + data[('megan', 'species', correction, simset)][0],
                  data[('megan', 'genus', correction, simset)][1] + data[('megan', 'genus', correction, simset)][0],
                  data[('megan', 'family', correction, simset)][1] + data[('megan', 'family', correction, simset)][0])  # false
        tp_megan = (data[('megan', 'species', correction, simset)][0], data[('megan', 'genus', correction, simset)][0],
                    data[('megan', 'family', correction, simset)][0])

        # Taxator-tk (species, genus, family)
        u_taxator = (100.0, 100.0, 100.0)
        w_taxator = (data[('taxator', 'species', correction, simset)][1] + data[('taxator', 'species', correction, simset)][0],
                    data[('taxator', 'genus', correction, simset)][1] + data[('taxator', 'genus', correction, simset)][0],
                    data[('taxator', 'family', correction, simset)][1] + data[('taxator', 'family', correction, simset)][0])  # false
        tp_taxator = (data[('taxator', 'species', correction, simset)][0], data[('taxator', 'genus', correction, simset)][0],
                      data[('taxator', 'family', correction, simset)][0])

        # Kraken (species, genus, family)
        u_kraken = (100.0, 100.0, 100.0)
        w_kraken = (data[('kraken', 'species', correction, simset)][1] + data[('kraken', 'species', correction, simset)][0],
                    data[('kraken', 'genus', correction, simset)][1] + data[('kraken', 'genus', correction, simset)][0],
                    data[('kraken', 'family', correction, simset)][1] + data[('kraken', 'family', correction, simset)][0])  # false
        tp_kraken = (data[('kraken', 'species', correction, simset)][0], data[('kraken', 'genus', correction, simset)][0],
                      data[('kraken', 'family', correction, simset)][0])

        # Taxator
        p1 = ax.bar(ind, u_taxator, width, color='y')
        p2 = ax.bar(ind, w_taxator, width, color='r')
        p3 = ax.bar(ind, tp_taxator, width, color='b')

        print('Taxator: %s %s %s' % (simset, correction, getMeasures(u_taxator, w_taxator, tp_taxator)))

        # PPS
        p1 = ax.bar(ind + width2, u_pps, width, color='y')
        p2 = ax.bar(ind + width2, w_pps, width, color='r')
        p3 = ax.bar(ind + width2, tp_pps,   width, color='b')

        print('PPS: %s %s %s' % (simset, correction, getMeasures(u_pps, w_pps, tp_pps)))

        # Megan
        p1 = ax.bar(ind + 2 * width2, u_megan, width, color='y')
        p2 = ax.bar(ind + 2 * width2, w_megan, width, color='r')
        p3 = ax.bar(ind + 2 * width2, tp_megan,   width, color='b')

        print('Megan: %s %s %s' % (simset, correction, getMeasures(u_megan, w_megan, tp_megan)))

        # Kraken
        p1 = ax.bar(ind + 3 * width2, u_kraken, width, color='y')
        p2 = ax.bar(ind + 3 * width2, w_kraken, width, color='r')
        p3 = ax.bar(ind + 3 * width2, tp_kraken, width, color='b')

        print('Kraken: %s %s %s' % (simset, correction, getMeasures(u_kraken, w_kraken, tp_kraken)))

        # PPS+
        p1 = ax.bar(ind + 4 * width2, u_ppsp, width, color='y')
        p2 = ax.bar(ind + 4 * width2, w_ppsp, width, color='r')
        p3 = ax.bar(ind + 4 * width2, tp_ppsp, width, color='b')

        print('PPS+: %s %s %s' % (simset, correction, getMeasures(u_ppsp, w_ppsp, tp_ppsp)))

        ax.set_ylabel('Percentage')
        ax.set_yticks(np.arange(0, 110, 10))

        ax.set_title(title + '\n')

        ticks = (width/2, width2 + width/2, width2*2 + width/2, width2*3 + width/2, width2*4 + width/2)
        ax.set_xticks(ticks + tuple(map(lambda x: x + 1, ticks)) + tuple(map(lambda x: x + 2, ticks)))
        ax.set_xticklabels(('taxator-tk', 'PPS', 'MEGAN', 'Kraken', 'PPS+',
                            'taxator-tk', 'PPS', 'MEGAN', 'Kraken', 'PPS+',
                            'taxator-tk', 'PPS', 'MEGAN', 'Kraken', 'PPS+'), rotation='vertical')

        ax.annotate(letter, xy=(1, 2), xytext=(-0.3, 115.))

        #ax.set_xticks((0.5,1.5,2.5), minor=True)
        #ax.set_xticklabels(['species', 'genus', 'family'], minor=True)  # how to put this up?

        #ax.xaxis.tick_top()

        ax.set_xlabel('species      genus        family  ')
        ax.xaxis.set_label_position('top')

    # add legend
    fontP = FontProperties()
    fontP.set_size(12)
    plt.legend((p3[0], p2[0], p1[0]), ('Correct', 'Incorrect', 'Unassigned'), loc='lower right',
                bbox_to_anchor=(0.57, -0.77), prop=fontP, ncol=3)

    plt.savefig(outFile, dpi=dpi)


def getMeasures(uList, wList, tpList):

    ret = []
    for u, w, tp in zip(uList, wList, tpList):

        correct = tp
        false = w - tp
        unassigned = u - w

        if not np.isclose(correct + false, 0.):
            p = round((correct / (correct + false)) * 100.0, 1)  # precision
        else:
            p = 0.0

        if not np.isclose(correct + false + unassigned, 0.):
            r = round((correct / (correct + false + unassigned)) * 100.0, 1)  # recall
        else:
            r = 0.0

        if not np.isclose(p + r, 0.):
            f = round(((2 * p * r) / (p + r)), 1)  # f-score, harmonic mean: (2 * precision * recall) / (precision + recall)
        else:
            f = 0.0
        ret.append('f:%s, p:%s, r:%s, (%s, %s, %s)' % (f, p, r, round(correct, 1), round(false, 1), round(unassigned, 1)))

    return ret


def createDemo():
    """
        Displays a demo chart containing one case (out of 4).
    """
    N = 3

    # PPS+ (species, genus, family)
    u_ppsp = (100.0, 100.0, 100.0)
    w_ppsp = (80.0, 60.0, 65.0)
    tp_ppsp = (20.0, 40.0, 50.0)

    # PPS (species, genus, family)
    u_pps = (100.0, 100.0, 100.0)
    w_pps = (80.0, 60.0, 65.0)
    tp_pps = (20.0, 40.0, 50.0)

    # Megan (species, genus, family)
    u_megan = (100.0, 100.0, 100.0)
    w_megan = (80.0, 60.0, 65.0)
    tp_megan = (20.0, 40.0, 50.0)

    # Taxator-tk (species, genus, family)
    u_taxator = (100.0, 100.0, 100.0)
    w_taxator = (80.0, 60.0, 65.0)
    tp_taxator = (20.0, 40.0, 50.0)

    ind = np.arange(N)    # the x locations for the groups
    width = 0.2
    width2 = 0.22

    # PPS+
    p1 = plt.bar(ind, u_ppsp, width, color='y')
    p2 = plt.bar(ind, w_ppsp, width, color='r')
    p3 = plt.bar(ind, tp_ppsp, width, color='b')

    # PPS
    p1 = plt.bar(ind + width2, u_pps, width, color='y')
    p2 = plt.bar(ind + width2, w_pps, width, color='r')
    p3 = plt.bar(ind + width2, tp_pps,   width, color='b')

    # Megan
    p1 = plt.bar(ind + 2 * width2, u_megan, width, color='y')
    p2 = plt.bar(ind + 2 * width2, w_megan, width, color='r')
    p3 = plt.bar(ind + 2 * width2, tp_megan,   width, color='b')

    # Taxator
    p1 = plt.bar(ind + 3 * width2, u_taxator, width, color='y')
    p2 = plt.bar(ind + 3 * width2, w_taxator, width, color='r')
    p3 = plt.bar(ind + 3 * width2, tp_taxator, width, color='b')

    plt.ylabel('Percentages')
    plt.title('Overall assignments \nspecies                           genus                             family      ')

    ticks = (width/2, width2 + width/2, width2*2 + width/2, width2*3 + width/2)
    plt.xticks(ticks + tuple(map(lambda x: x + 1, ticks)) + tuple(map(lambda x: x + 2, ticks)),
    ('PPS+', 'PPS', 'Megan', 'taxator-tk', 'PPS+', 'PPS', 'Megan', 'taxator-tk', 'PPS+', 'PPS', 'Megan', 'taxator-tk',), rotation=90)
    plt.yticks(np.arange(0, 110, 10))
    plt.legend( (p1[0], p2[0], p3[0]), ('Unassigned', 'False', 'True'), loc='upper left' )
    plt.subplots_adjust(left=0.1, right=0.9, top=0.90, bottom=0.15)

    plt.show()

def assignedAtRank(rankR, mapFile, fasType):
    if mapFile is None or not os.path.isfile(mapFile):
        print('Nothing %s%s%s' % (rankR, mapFile, fasType))
        return 0

    # mapFile = '/Users/ivan/Documents/work/Doc/PPSplus/revision/kraken/results/kraken_lognorm_genus_lab.csv'
    if fasType == 'lognorm':
        fasFile = '/Users/ivan/Documents/work/Doc/PPSplus/revision/nobackup/VelvetNoscafMinimus2_ma_merge_lognorm_min_1000bp.fna'
    else:
        assert fasType == 'uniform'
        fasFile = '/Users/ivan/Documents/work/Doc/PPSplus/revision/nobackup/VelvetNoscafMinimus2_ma_merge_uniform_min_1000bp.fna'

    taxFile = '/Users/ivan/Documents/work/Doc/PPSplus/revision/nobackup/ncbitax_sqlite.db'

    predDict = csv.predToDict(mapFile)
    bpDict = fas.getSequenceToBpDict(fasFile)
    tax = taxonomy_ncbi.TaxonomyNcbi(taxFile)

    rankList = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
    rankSet = set(rankList)
    rankMap = {}
    for i in range(len(rankList)):
        rankMap[rankList[i]] = i

    counts = np.zeros(len(rankList), np.int64)

    for seqId, taxonId in predDict.iteritems():
        bp = bpDict[seqId]

        rank = tax.getRank(taxonId)
        while rank not in rankSet and taxonId != 1 and taxonId is not None:
            taxonId = tax.getParentNcbid(taxonId)
            rank = tax.getRank(taxonId)
            # print taxonId

        if rank in rankSet:
            counts[rankMap[rank]] += bp

    # print counts

    for i in range(len(rankList) -1, 0, -1):
        counts[i-1] += counts[i]

    # print counts

    # print ','.join(map(lambda x: '%s:%s' % (x[0], x[1]), zip(rankList, counts)))
    d = dict(map(lambda x: (x[0], x[1]), zip(rankList, counts)))

    tax.close()

    return d[rankR]

if __name__ == "__main__":
    # print('DO NOT remove all nodes of degree 2, only if the record is redundant!!!')
    createChart()
    #createDemo()

