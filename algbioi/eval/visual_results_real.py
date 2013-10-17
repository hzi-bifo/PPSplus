"""
    To visualize results for real datasets.

    contigs vs scaffolds: /Users/ivan/Documents/work/Doc/PPSplus/evaluation_real/contigs_vs_scaff


"""
import os
import math

import numpy as np
import matplotlib.pyplot as plt

# settings
CONTIGS_VS_SCAFFOLDS_INPUT_DIR = '/Users/ivan/Documents/work/Doc/PPSplus/evaluation_real/contigs_vs_scaff'
CONTIGS_VS_SCAFFOLDS_OUTPUT_DIR = '/Users/ivan/Documents/work/Doc/PPSplus/figures_real/contigs_vs_scaff'


def heatMapPlot(dataModel, outFile=None, dpi=300, showFigure=False, kbp=True, plotUnassigned=True):
    """
        Plots a heat map.

        @param dataModel: contains data to plot the chart
        @type dataModel: HeatMapDataModel
        @param outFile: output file to store the chart
        @param dpi: output file dpi
        @param showFigure: show the figure if true
        @param kbp: plot data in terms of bp (else in terms of count)
    """
    # data in terms of bp or count
    if kbp:
        data = dataModel.getMatrixKBp()
    else:
        data = dataModel.getMatrixCount()

    # headers
    rowHeader = dataModel.getRowHeader()
    colHeader = dataModel.getColHeader()

    if not plotUnassigned:
        tmp = []
        for row in data[:-1]:
            tmp.append(row[:-1])
        data = tmp
        rowHeader = rowHeader[:-1]
        colHeader = colHeader[:-1]

    if len(data) == 0:
        print("No data to be plotted for: %s" % outFile)
        return

    # normalize data
    max = 0
    for row in data:
        for col in row:
            if col > max:
                max = col

    base = max + 1
    nF = lambda x: math.log(1 + x, base)
    normData = []

    for row in data:
        r = []
        for col in row:
            r.append(nF(col))
        normData.append(r)
    data = np.array(normData)

    # row and column headers
    if len(rowHeader) > 40:
        rowHeader = len(rowHeader) * ['']
    if len(colHeader) > 40:
        colHeader = len(colHeader) * ['']

    # create heat map
    fig, ax = plt.subplots()
    ax.pcolor(data, cmap=plt.cm.Reds)

    # set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=5)
    plt.tick_params(axis='both', which='minor', labelsize=5)

    # set figure size depending on the amount of the data
    height = len(rowHeader) * 0.15
    width = len(colHeader) * 0.15

    # assure min size

    if height < 3.:
        height = 3.
    if width < 3.:
        width = 3.

    if height > 6.:
        height = 6.
    if width > 6.:
        width = 6.

    fig.set_size_inches(width, height)

    # adjust spacing
    left = 0.125
    right = 0.9
    bottom = 0.1
    top = 0.9
    wspace = 0.2
    hspace = 0.2
    #xlabelFontSize = 10

    if len(rowHeader) < 5:
        top = 0.8
        bottom = 0.1

    elif len(rowHeader) < 15:
        top = 0.6
        bottom = 0.1

    elif len(rowHeader) < 20:
        top = 0.65
        bottom = 0.1

    elif len(rowHeader) < 35:
        top = 0.69
        bottom = 0.1
    else:
        top = 0.8
        bottom = 0.1
####
    if len(colHeader) < 5:
        right = 0.9
        left = 0.2
        xlabelFontSize = 4

    elif len(colHeader) < 15:
        right = 0.85
        left = 0.4
        xlabelFontSize = 4

    elif len(colHeader) < 20:
        right = 0.9
        left = 0.3
        xlabelFontSize = 4

    elif len(colHeader) < 25:
        right = 0.9
        left = 0.3
        xlabelFontSize = 5

    elif len(colHeader) < 35:
        right = 0.9
        left = 0.3
        xlabelFontSize = 6


    elif len(colHeader) < 40:
        right = 0.95
        left = 0.2
        xlabelFontSize = 7

    elif len(colHeader) < 60:
        right = 0.95
        left = 0.3
        xlabelFontSize = 9
    else:
        right = 0.95
        left = 0.3
        xlabelFontSize = 9


    if len(rowHeader) > 40 and len(colHeader) > 40:
        top = 0.9
        left = 0.1
    elif len(colHeader) > 40:
        top = 0.9

    plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

    #plt.subplots_adjust(top=0.8, bottom=0.12, right=0.95)


    # put the major ticks at the middle of each cell, notice "reverse" use of dimension
    ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

    # invert the y axis, x ticks are on top and rotated
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    plt.xticks(rotation=90)

    # axis limits
    ax.set_xlim(-0.1, len(colHeader) + 0.1)
    ax.set_ylim(len(rowHeader) + 0.1, -0.1)

    # set axis labels
    ax.set_xticklabels(colHeader, minor=False)
    ax.set_yticklabels(rowHeader, minor=False)

    # set chart text - info
    if dataModel.getInfo() is not None:
        #ax.set_xlabel(dataModel.getInfo())
        plt.xlabel(dataModel.getInfo(), fontsize=xlabelFontSize)
        #ax.text(0, len(rowHeader) + len(colHeader) / 30., dataModel.getInfo(), fontsize=len(colHeader) / 4.)  # len(colHeader) / 3, len(rowHeader) + 4
        # style='italic' bbox={'facecolor':'white', 'alpha':0.5, 'pad':10}

    # save figure to file
    if outFile is not None:
        plt.savefig(outFile, dpi=dpi)

    # show figure
    if showFigure:
        plt.show()

    # clear the current plot (enables to start a new plot)
    plt.clf()


class HeatMapDataModel:
    def __init__(self):
        """
            Holds data model for a heat map.
        """
        self.__rowHeader = None
        self.__colHeader = None
        self.__matrixKBp = None
        self.__matrixCount = None
        self.__info = None

    def rowHeader(self, rowHeader):
        self.__rowHeader = rowHeader
        return self

    def colHeader(self, colHeader):
        self.__colHeader = colHeader
        return self

    def matrixKBp(self, matrix):
        self.__matrixKBp = matrix
        return self

    def matrixCount(self, matrix):
        self.__matrixCount = matrix
        return self

    def info(self, info):
        self.__info = info
        return self

    def getColHeader(self):
        return self.__colHeader

    def getRowHeader(self):
        return self.__rowHeader

    def getMatrixKBp(self):
        return self.__matrixKBp

    def getMatrixCount(self):
        return self.__matrixCount

    def getInfo(self):
        return self.__info


class CTDataHandler:
    def __init__(self, inputDir):
        """
            Holds all data contained in the input directory where the input directory contains subdirectories with the
            csv files containing comparison tables.
        """
        self.__subDictToEntryList = {}
        for d in os.listdir(inputDir):
            dPath = os.path.join(inputDir, d)
            l = []
            for f in os.listdir(dPath):
                l.append(CTEntryHandler(os.path.join(dPath, f)))
                #break  # test
            self.__subDictToEntryList[d] = l
            #break  # test

    def generator(self):
        """
            Returns pairs: (dirName, list of entries)
            @rtype: (str, list of CTEntryHandler)
        """
        for k, v in self.__subDictToEntryList.iteritems():
            yield (k, v)


class CTEntryHandler:
    def __init__(self, entryFilePath):
        """
            Represents one comparison table / heat map.

            @param entryFilePath:
        """
        self.__dataModel = HeatMapDataModel()

        # read in the data from the file as a matrix
        matrix = []
        for line in open(entryFilePath):
            matrix.append(map(lambda x: x.strip(), line.split(',')))

        # set header (first line)
        self.__dataModel.colHeader(matrix[0][1:])

        # set row header (first column until empty line)
        rowHeader = []
        for row in matrix[1:]:
            if row[0] == '':
                break
            else:
                rowHeader.append(row[0])
        self.__dataModel.rowHeader(rowHeader)

        # read in the entries in the matrix (as kbp and as count)
        matrixKBp = []
        matrixCount = []

        for row in matrix[1:]:
            if row[0] == '':
                break
            rowKBp = []
            rowCount = []
            for col in row[1:]:
                kbp = 0
                count = 0
                if col != '':
                    kbp = int(col.split('k')[0])
                    count = int(col.split('(')[1].rstrip(')'))
                rowKBp.append(kbp)
                rowCount.append(count)
            matrixKBp.append(rowKBp)
            matrixCount.append(rowCount)

        self.__dataModel.matrixKBp(matrixKBp).matrixCount(matrixCount)

        # read matched entries
        meta = {}
        for row in matrix:
            if 'Matches' in row[0]:
                meta['Matches'] = map(lambda x: float(x.strip('k%')), row[1:])
                break

        # generate heat map info
        dataset, method = os.path.dirname(entryFilePath).split(os.sep)[-1].split('_')[0:2]
        rank = os.path.basename(entryFilePath).split('.')[-2].split('_')[0]
        methodMap = { 'megan': 'Megan', 'taxator': 'Taxator-tk', 'pps': 'PhyloPythiaS', 'ppsp': 'PhyloPythiaS plus'}
        datasetMap = {'cr': 'Cow Rumen', 'hg': 'Human Gut'}

        info = 'Consistency, %s, %s, %s, agree: %skb %s%% ' % \
               (methodMap[method], datasetMap[dataset], rank, intPrint(int(meta['Matches'][0])), meta['Matches'][2])

        self.__dataModel.info(info)
        self.__outFilePrefix = '%s_%s_%s' % (dataset, method, rank)

    def plotEntry(self, outDir, showFigure=False, dpi=300, plotUnassigned=True):
        heatMapPlot(self.__dataModel, showFigure=showFigure, dpi=dpi,
                    outFile=os.path.join(outDir, str(self.__outFilePrefix) + '.png'), plotUnassigned=plotUnassigned)


def intPrint(i):
    """
        @type i: int
        @rtype: str
    """
    a = list(str(i))
    a.reverse()
    ret = ''
    k = 0
    for e in a:
        ret = e + ret
        k += 1
        if k % 3 == 0:
            ret = ',' + ret
    return ret.strip(',')


def _main():
    for directory, entryList in CTDataHandler(CONTIGS_VS_SCAFFOLDS_INPUT_DIR).generator():
        outDir = os.path.join(CONTIGS_VS_SCAFFOLDS_OUTPUT_DIR, directory)
        if not os.path.isdir(outDir):
            os.mkdir(outDir)
        #if 'taxator' not in outDir:
        #    continue
        for entry in entryList:
            #try:
            entry.plotEntry(outDir, plotUnassigned=False, dpi=300)
            #except Exception as e:
            #    print('exception: ' + e.message + ' ' + outDir)



def _test():
    for k, v in CTDataHandler(CONTIGS_VS_SCAFFOLDS_INPUT_DIR).generator():
        v[0].plotEntry()
        break


if __name__ == "__main__":
    #_test()
    _main()