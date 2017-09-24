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

    Miscellaneous functions to test the SILVA DATABASE
"""

from algbioi.com import taxonomy_ncbi
from algbioi.com import csv


class TaxonomyWrapper():
    def __init__(self, taxonomy):
        """
            @type taxonomy: taxonomy_ncbi.TaxonomyNcbi
        """
        self._toParent = {}
        self._toRank = {}
        self._taxonomy = taxonomy

    def getParentNcbid(self, taxid):
        if taxid in self._toParent:
            return self._toParent[taxid]
        else:
            parent = self._taxonomy.getParentNcbid(taxid)
            self._toParent[taxid] = parent
            return parent

    def getRank(self, taxid):
        if taxid in self._toRank:
            return self._toRank[taxid]
        else:
            rank = self._taxonomy.getRank(taxid)
            self._toRank[taxid] = rank
            return rank

    def close(self):
        self._taxonomy.close()


def test01():
    peter = '/Volumes/VerbatimSSD/work/vm_1_4/peter/SSURef_NR99_115.filter_bac_arc.map.csv'
    me = '/Volumes/VerbatimSSD/work/vm_1_4/reference_NCBI20140513/silva115/db/ssu_115.tax'
    db = '/Volumes/VerbatimSSD/work/vm_1_4/reference_NCBI20140513/taxonomy/ncbitax_sqlite.db'
    t = TaxonomyWrapper(taxonomy_ncbi.TaxonomyNcbi(db))

    p = csv.predToDict(peter)

    m = {}
    ml = {}
    for line in open(me):
        k, v = line.strip().strip(';').split('\t')

        k = k.split('|')[0].replace('-', '.')
        ml[k] = v
        m[k] = set(map(lambda x: int(x), v.split(';')))

    # sameNameCount = 0
    labelConsistent = 0
    labelNotConsistent = 0
    # counter = 0
    dRank = {}
    for k, v in p.iteritems():
        vv = v
        # counter += 1
        # if counter % 1000 == 0:
        #     print counter
        if k in m:
            # sameNameCount += 1
            if v in m[k]:
                labelConsistent += 1
            else:
                rank = t.getRank(v)
                if rank == 'species':
                    labelNotConsistent += 1
                    print k, ml[k], vv
                else:
                    while v != 1 and rank != 'species' and v is not None:
                        v = t.getParentNcbid(v)
                        rank = t.getParentNcbid(v)
                        # print rank
                    if rank == 'species' and v not in m[k]:
                        labelNotConsistent += 1
                        print k, m[k], vv
                    else:
                        labelConsistent += 1
        else:
            # print k, v
            rank = t.getRank(v)
            if rank in dRank:
                dRank[rank] += 1
            else:
                dRank[rank] = 0

    # print v, m[k]

    #     if k not in p:
    #         print '.'
    print 'len(m):', len(m)
    print 'len(p):', len(p)
    print 'labelConsistent:', labelConsistent
    print 'labelNotConsistent', labelNotConsistent

    for k, v in dRank.iteritems():
        print k, v

    t.close()
    pass


if __name__ == "__main__":
    test01()