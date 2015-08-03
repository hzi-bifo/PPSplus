
#!/usr/bin/env python

"""
    Copyright (C) 2015  Ivan Gregor

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

    File preprocession for the container integration.
"""


import os
from algbioi.com import csv
import sys


def _main():
    if len(sys.argv) < 3:
        print('Parameters:\npipelineDir\ndstBinFile\n')

    pipelineDir, dstBinFile = sys.argv[1:3]

    assert os.path.isdir(pipelineDir)
    assert os.path.isdir(os.path.dirname(dstBinFile))

    workingDir = os.path.join(pipelineDir, 'working')
    for f in os.listdir(workingDir):
        fPath = os.path.join(workingDir, f)
        if os.path.isfile(fPath) and f.endswith('.binning'):
            out = csv.OutFileBuffer(dstBinFile)
            for line in open(fPath):
                line = line.strip()
                out.writeText('%s\n' % line)
            out.close()
            break

if __name__ == "__main__":
    _main()