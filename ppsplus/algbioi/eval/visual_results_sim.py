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


    Provide methods to visualize evaluation results.
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from algbioi.com import taxonomy_ncbi


def simSetPlot(dataModel, outFile=None, dpi=300, showFigure=True, taxonomicRanks=taxonomy_ncbi.TAXONOMIC_RANKS[1:]):
    """
        Plot precision and recall for simulated datasets.
        
        @param dataModel: data to plot
        @type dataModel: SimSetDataModel
        @param outFile: output file to store results
        @param dpi: output file dpi
        @param showFigure: show the figure if true
    """
    # set plot title
    plt.title(dataModel.getTitle())

    # set x axis values
    x = range(1, len(taxonomicRanks) + 1)

    # set line properties
    lineWidth = 1.5
    marker = 'o'
    markerSize = 3
    colors = ['r', 'b', 'g', 'c', 'm']  # red, blue, green, cyan, magenta
    lineStylePrecision = '-'  # solid line style
    lineStyleRecall = '--'  # dashed line style

    # precision lines
    for i in range(len(dataModel.getPrecisionList())):
        precisionList = dataModel.getPrecisionList()[i]
        if len(precisionList) < len(taxonomicRanks):
            precisionList += (len(taxonomicRanks) - len(precisionList)) * [0.0]
        plt.plot(x, precisionList, linestyle=lineStylePrecision, linewidth=lineWidth, color=colors[i], marker=marker,
                 markersize=markerSize, label='(P) ' + dataModel.getLegend()[i])

    # recall lines
    for i in range(len(dataModel.getRecallList())):
        recallList = dataModel.getRecallList()[i]
        if len(recallList) < len(taxonomicRanks):
            recallList += (len(taxonomicRanks) - len(recallList)) * [0.0]
        plt.plot(x, recallList, linestyle=lineStyleRecall, linewidth=lineWidth, color=colors[i], marker=marker,
                 markersize=markerSize, label='(R) ' + dataModel.getLegend()[i])

    # set legend (http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot)
    legend = plt.legend(title=dataModel.getLegendTitle(), loc='lower left', fontsize='medium', ncol=2)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')  # set legend color

    # x axis ticks (You can specify a rotation for the tick labels in degrees or with keywords.)
    plt.xticks(x, taxonomicRanks, rotation=20)

    # margin (Pad margins so that markers don't get clipped by the axes)
    plt.margins(0.2)

    # tick labels spacing (Tweak spacing to prevent clipping of tick-labels)
    plt.subplots_adjust(bottom=0.15)

    # grid
    ax = plt.axes()
    ax.xaxis.grid()  # vertical lines
    ax.yaxis.grid()  # horizontal lines

    # axis limits
    ax.set_xlim(0.9, 7.1)
    ax.set_ylim(0, 1.02)

    # y-axis label formatting
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos=0: '%1.f%%' % (100*y)))

    # y-axis ticks frequency
    plt.yticks(np.arange(0.0, 1.1, 0.1))

    # save figure to a file
    if outFile is not None:
        plt.savefig(outFile, dpi=dpi)

    # show figure
    if showFigure:
        plt.show()

    # clear the current plot (enables to start a new plot)
    plt.clf()


class SimSetDataModel():
    def __init__(self):
        """
            Represents a data model for a simulated data set.
            (Initialize using the builder pattern.)
        """
        self.__title = None
        self.__precisionList = None
        self.__recallList = None
        self.__legend = None
        self.__legendTitle = None

    def title(self, title):
        """
            @param title: Figure title
        """
        self.__title = title
        return self

    def precision(self, precisionList):
        """
            @param precisionList: two dimensional array, each row represents precision vector from the top level rank
        """
        self.__precisionList = precisionList
        return self

    def recall(self, recallList):
        """
            @param recallList: two dimensional array, each row represents recall vector from the top level rank
        """
        self.__recallList = recallList
        return self

    def legend(self, legend):
        """
            @param legend: a list of legends one for each line, first precision, then recall
        """
        self.__legend = legend
        return self

    def legendTitle(self, legendTitle):
        """
            @param legendTitle: title of the legend box
        """
        self.__legendTitle = legendTitle
        return self

    def getTitle(self):
        return self.__title

    def getPrecisionList(self):
        return self.__precisionList

    def getRecallList(self):
        return self.__recallList

    def getLegend(self):
        return self.__legend

    def getLegendTitle(self):
        return self.__legendTitle


class SimFileHandler():
    def __init__(self, simFile):
        """
            Holds precision/recall information for one Method/basic settings.
            Ignore lines starting with %.
            Lines starting with # are considered as headers
            Empty lines are considered as record separators.

            @param simFile: input file with all the precision/recall information
        """
        self.__results = []
        handler = SimFileEntryHandler()

        for line in open(simFile):
            if line.startswith('%'):
                continue
            if line.strip() != "":
                handler.parseLine(line)
            else:
                if not handler.isEmpty():
                    self.__results.append(handler)
                handler = SimFileEntryHandler()

    def getResultList(self):
        """
            Gets all entries as a list result.
        """
        return self.__results

    def generator(self):
        """
            Gets all entries as a generator.
            @rtype: list of SimFileEntryHandler
        """
        for result in self.__results:
            yield result


class SimFileEntryHandler():
    def __init__(self):
        """
            Holds information about one entry.
        """
        self.__headers = ""
        self.__ranksSet = set(taxonomy_ncbi.TAXONOMIC_RANKS)
        self.__weightedBpDict = {}
        self.__notWeightedBpDict = {}
        self.__rank = None

    def parseLine(self, line):
        if line.startswith('#'):
            self.__headers += line
            #self.__rank = None
        elif line.split(',')[0] in self.__ranksSet:
            self.__rank = line.split(',')[0]
        elif self.__rank is not None:
            lineTokens = line.split(',')

            if 'bp' in lineTokens[4]:
                if 'weighted' in lineTokens[5]:
                    self.__weightedBpDict[self.__rank] = line.split(',')
                elif 'not weighted' in lineTokens[5]:
                    self.__notWeightedBpDict[self.__rank] = line.split(',')
                else:
                    print('SimResultHandler.parseLine: wrong token')

    def getEntryList(self, field='precision', weighted=True, bp=True):
        """
            To get a list of precision or recall values.

            @param field: precision or recall
            @type field: str
            @param weighted: weighted variant of the precision and recall values
            @param bp: must be true
            @return: a list starting with a value for the highest rank
        """
        assert bp
        if 'precision' in field:
            idx = 0
        elif 'recall' in field:
            idx = 1
        else:
            raise ValueError('The parameter value can be only precision or recall!')
        retList = []
        for rank in taxonomy_ncbi.TAXONOMIC_RANKS:
            if weighted:
                if rank in self.__weightedBpDict:
                    retList.append(self.__weightedBpDict[rank][idx])
            else:
                if rank in self.__notWeightedBpDict:
                    retList.append(self.__notWeightedBpDict[rank][idx])
        return retList

    def headerContains(self, keywordList):
        """
            Returns true if the header contains all the entries from the keywordList.

            @param keywordList: a list of keywords
            @type keywordList: list of str
            @rtype: bool
        """
        for keyword in keywordList:
            if keyword not in self.__headers:
                return False
        return True

    def isEmpty(self):
        if (len(self.__weightedBpDict) == 0) and (len(self.__notWeightedBpDict) == 0):
            return True
        else:
            return False


class SimPlotSettings():
    def __init__(self, fileDir, fileKeywords, scenarioList, title, legend, legendTitle):
        """
            Represents one chart for the simulated data sets scenario.
            One chart can represent results for several settings in one plot.
            The results file containing given keywords will be taken as an input file.
            An entry in the results file will be taken if its header contain given keywords.

            @param fileDir: a directory containing results files
            @param fileKeywords: a list of keywords that a line of a file starting with % must contain to be identified
            @param scenarioList: a list of lists where each list represents a scenario
            @param title: chart title
            @param legend: entries of the legend (will be used twice, for precision and recall)
            @param legendTitle: title of the whole legend
            @param fileDir: output chart file
        """
        # find the right file
        inFile = None
        for f in os.listdir(fileDir):
            fPath = os.path.join(fileDir, f)
            for line in open(fPath):
                if line.startswith('%'):
                    ok = True
                    for keyword in fileKeywords:
                        if keyword not in line:
                            ok = False
                    if ok:
                        inFile = fPath
                if inFile is not None:
                    break
            if inFile is not None:
                break
        if inFile is None:
            raise ValueError("Can't find file identified by %s %s" % (fileDir, fileKeywords))

        # Read in data
        inFileHandler = SimFileHandler(inFile)
        precisionList = []
        recallList = []
        toFloat = lambda x: float(x) / 100.  # to map percentages to fractions

        # For each scenario, get precision and recall values
        for scenarioKeywords in scenarioList:
            for handler in inFileHandler.generator():
                if handler.headerContains(scenarioKeywords):
                    precisionList.append(map(toFloat, handler.getEntryList(field='precision', weighted=True)))
                    recallList.append(map(toFloat, handler.getEntryList(field='recall', weighted=True)))

        # Builds data model.
        self.__dataModel = SimSetDataModel().precision(precisionList).recall(recallList).\
            title(title).legend(legend).legendTitle(legendTitle)

    def generateChart(self, outFile, showFigure=False):
        """
            Generates the chart.

            @param outFile: output file to store the chart
            @param showFigure: show the figure if true
        """
        simSetPlot(self.__dataModel, outFile, showFigure=showFigure)


def melanieSimsetGenerateCharts():
    """
        This method generates charts for the simulated datasets from Melanie.
        Datasets were generated with two distributions: uniform and log-normal.
        Binning methods used are: PPSP, PPS, Megan, Taxator

    """
    # in/out directories
    # fileDir = '/Users/ivan/Documents/work/Doc/PPSplus/evaluation_simset'
    fileDir = '/Users/ivan/Documents/work/Doc/PPSplus/revision/evaluation_simset'
    # outFig = '/Users/ivan/Documents/work/Doc/PPSplus/figures_simset'
    outFig = '/Users/ivan/Documents/work/Doc/PPSplus/revision/figures_simset'

    # For Clark and Kraken
    for fileKeywords, simsetType, method, rank in zip(
            [['UNIFORM', 'kraken'], ['UNIFORM', 'clark'],
                ['LOG-NORM', 'kraken'], ['LOG-NORM', 'clark']],
            (2 * ['uniform'] + 2 * ['log-norm']),
            ['Kraken', 'CLARK', 'Kraken', 'CLARK'],
            ['species', 'species', 'species', 'species']):

        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['no exclude', 'nc'],
                        ['new strain', 'nc'],
                        ['new species', 'nc'],
                        ['new genus', 'nc']
                        ],
                        title='%s (%s dataset)' % (method, simsetType),
                        legend=['known strain', 'new strain', 'new species', 'new genus'],
                        legendTitle='precision (P), recall (R)').\
            generateChart(os.path.join(outFig, '%s_simset_%s_%s_rs_nc.png' % \
                                               (str.lower(method), simsetType, rank)))

        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['no exclude', 'correction'],
                        ['new strain', 'correction'],
                        ['new species', 'correction'],
                        ['new genus', 'correction']
                        ],
                        title='%s (%s dataset, correction)' % (method, simsetType),
                        legend=['known strain', 'new strain', 'new species', 'new genus'],
                        legendTitle='precision (P), recall (R)').\
            generateChart(os.path.join(outFig, '%s_simset_%s_%s_rs_correction.png' % \
                                               (str.lower(method), simsetType, rank)))

    # return  # TODO: !!! run only for Clark and kraken


    # For PPSP (modeled down to species)
    for fileKeywords, simsetType in zip([['UNIFORM', 'ppsp-species'], ['LOG-NORM', 'ppsp-species']],
                                       ['uniform', 'log-norm']):
        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['rs_no_mg_no', 'nc'],
                        ['rs_strain_mg_no', 'nc'],
                        ['rs_species_mg_no', 'nc'],
                        ['rs_genus_mg_no', 'nc']
                        ],
                        title='PhyloPythiaS plus (%s dataset)' % simsetType,
                        legend=['known strain', 'new strain (rs)', 'new species (rs)', 'new genus (rs)'],
                        legendTitle='Known strain (mg), precision (P), recall (R)').\
            generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_no_nc.png' % simsetType))

        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                    scenarioList=[
                    ['rs_no_mg_no', 'correction'],
                    ['rs_strain_mg_no', 'correction'],
                    ['rs_species_mg_no', 'correction'],
                    ['rs_genus_mg_no', 'correction']
                    ],
                    title='PhyloPythiaS plus (%s dataset, correction)' % simsetType,
                    legend=['known strain', 'new strain (rs)', 'new species (rs)', 'new genus (rs)'],
                    legendTitle='Known strain (mg), precision (P), recall (R)').\
        generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_no_correction.png' % simsetType))

        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['rs_strain_mg_strain', 'nc'],
                        ['rs_species_mg_strain', 'nc'],
                        ['rs_genus_mg_strain', 'nc']
                        ],
                        title='PhyloPythiaS plus (%s dataset)' % simsetType,
                        legend=['new strain (rs)', 'new species (rs)', 'new genus (rs)'],
                        legendTitle='New strain (mg), precision (P), recall (R)').\
            generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_strain_nc.png' % simsetType))


        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['rs_strain_mg_strain', 'nc'],
                        ['rs_species_mg_strain', 'nc'],
                        ['rs_genus_mg_strain', 'nc']
                        ],
                        title='PhyloPythiaS plus (%s dataset)' % simsetType,
                        legend=['new strain (rs)', 'new species (rs)', 'new genus (rs)'],
                        legendTitle='New strain (mg), precision (P), recall (R)').\
            generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_strain_nc.png' % simsetType))


        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['rs_strain_mg_strain', 'correction'],
                        ['rs_species_mg_strain', 'correction'],
                        ['rs_genus_mg_strain', 'correction']
                        ],
                        title='PhyloPythiaS plus (%s dataset, correction)' % simsetType,
                        legend=['new strain (rs)', 'new species (rs)', 'new genus (rs)'],
                        legendTitle='New strain (mg), precision (P), recall (R)').\
            generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_strain_correction.png' % simsetType))



        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['rs_no_mg_no', 'nc'],
                        ['rs_strain_mg_strain', 'nc'],
                        ['rs_species_mg_species', 'nc'],
                        ['rs_genus_mg_genus', 'nc']
                        ],
                        title='PhyloPythiaS plus (%s dataset)' % simsetType,
                        legend=['known strain', 'new strain', 'new species', 'new genus'],
                        legendTitle='precision (P), recall (R)').\
            generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_rs_nc.png' % simsetType))


        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['rs_no_mg_no', 'correction'],
                        ['rs_strain_mg_strain', 'correction'],
                        ['rs_species_mg_species', 'correction'],
                        ['rs_genus_mg_genus', 'correction']
                        ],
                        title='PhyloPythiaS plus (%s dataset, correction)' % simsetType,
                        legend=['known strain', 'new strain', 'new species', 'new genus'],
                        legendTitle='precision (P), recall (R)').\
            generateChart(os.path.join(outFig, 'ppsp_simset_%s_species_mg_rs_correction.png' % simsetType))



    # For PPS, Megan, Taxator
    for fileKeywords, simsetType, method, rank in zip(
            [['UNIFORM', 'megan'], ['UNIFORM', 'taxator'], ['UNIFORM', 'pps-genus'],
                ['LOG-NORM', 'megan'], ['LOG-NORM', 'taxator'], ['LOG-NORM', 'pps-genus']],
            (3 * ['uniform'] + 3 * ['log-norm']),
            ['Megan', 'Taxator-tk', 'PhyloPythiaS', 'Megan', 'Taxator-tk', 'PhyloPythiaS'],
            ['species', 'species', 'genus', 'species', 'species', 'genus']):

        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['no exclude', 'nc'],
                        ['new strain', 'nc'],
                        ['new species', 'nc'],
                        ['new genus', 'nc']
                        ],
                        title='%s (%s dataset)' % (method, simsetType),
                        legend=['known strain', 'new strain', 'new species', 'new genus'],
                        legendTitle='precision (P), recall (R)').\
            generateChart(os.path.join(outFig, '%s_simset_%s_%s_rs_nc.png' % \
                                               (str.lower(method), simsetType, rank)))

        SimPlotSettings(fileDir, fileKeywords=fileKeywords,
                        scenarioList=[
                        ['no exclude', 'correction'],
                        ['new strain', 'correction'],
                        ['new species', 'correction'],
                        ['new genus', 'correction']
                        ],
                        title='%s (%s dataset, correction)' % (method, simsetType),
                        legend=['known strain', 'new strain', 'new species', 'new genus'],
                        legendTitle='precision (P), recall (R)').\
            generateChart(os.path.join(outFig, '%s_simset_%s_%s_rs_correction.png' % \
                                               (str.lower(method), simsetType, rank)))


def _test():
    yPrecision = [
        [.98, .98, .95, .92, .9, .88, .85],
        [.95, .85, .75, .65, .55, .45, .35],
        [.8, .75, .7, .65, .6, .55, .5],
        [.9, .8, .7, .6, .5, .4, .3]]
    yRecall = [
        [.88, .88, .85, .82, .8, .78, .75],
        [.9, .8, .75, .71, .65, .6, .54],
        [.7, .78, .7, .65, .64, .61, .58],
        [.6, .7, .6, .6, .5, .5, .5]]

    dataModel = SimSetDataModel().\
        title('PhyloPythiaS plus').\
        precision(yPrecision).recall(yRecall).\
        legend(['known strain', 'new strain (rs)', 'new species (rs)', 'new genus (rs)']).\
        legendTitle('New species (mg), precision (P), recall (R)')

    simSetPlot(dataModel, outFile='/Users/ivan/Documents/work/Doc/PPSplus/revision/figures_simset/test.png')

def _main():
    melanieSimsetGenerateCharts()

if __name__ == "__main__":
    _main()
    # _test()