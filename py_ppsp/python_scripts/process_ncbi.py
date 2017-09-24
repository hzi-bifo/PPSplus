#!/usr/bin/env python
"""
Author: Jessika Fiedler (jessika.fiedler@hhu.de)
Copyright (c) Jessika Fiedler 2014

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

About:
process ncbi download by creating fasta files with appropriate names
supports gbk and gbff files
"""

import sys
import os
import gbk


class Process_NCBI:
    
    def __init__(self, n_dir, o_dir, e=False):
        assert os.path.isdir(n_dir) and os.path.isdir(o_dir)

        self.ncbi_dir = n_dir
        self.output_dir = o_dir
        self.exclude = e

    def run(self):
        files = self.get_files_to_process()
        if len(files) == 0:
            sys.stderr.write("No files to process..")
            return False
        self.process_files(files)
        return True

    def get_files_to_process(self):
        process_extensions = ["gbk", "gbff"]
        sys.stdout.write("\nPreparing ...")

        ncbi_dir_contents = os.listdir(self.ncbi_dir)
        assert len(ncbi_dir_contents) > 0

        files_to_process = {}   # mapping file -> directory
        for ndc in ncbi_dir_contents:
            if os.path.isdir(os.path.join(self.ncbi_dir, ndc)):
                for c in os.listdir(os.path.join(self.ncbi_dir, ndc)):
                    if c.split('.')[-1] in process_extensions:
                        path = "{folder}/{subfolder}/{file}".format(folder=self.ncbi_dir, subfolder=ndc, file=c)
                        files_to_process[path] = ndc
            if ndc.split('.') in process_extensions:
                files_to_process[os.path.join(self.ncbi_dir, ndc)] = self.ncbi_dir

        return files_to_process

    def process_files(self, files_to_process):
        assert len(files_to_process) > 0

        sys.stdout.write("done\nProcessing (this will take a while)...")
        # A=Ask, S=Skip, O=Overwrite, SA=Skip All, OA=Overwrite All
        action_file_exists = "A"
        valid_answers_file_exist = ["S", "O", "SA", "OA"]
        question_file_exist = "The file exists what should we do [S:Skip(default),O:Overwrite,SA:Skip All,OA:Overwrite All]? "

        processed_taxon = set()
        processed_dir = set()

        for f, d in files_to_process.iteritems():
            records = gbk.readFromGbkFile(f, listOfKeywords=['description', 'accession', 'gi', 'taxonId', 'seq'])
            for record in records:
                if self.exclude and ('phage' in record['description'] or 'plasmid' in record['description']):
                    continue
                fastafile = "{outdir}{sep}{taxid}.1.fna".format(outdir=self.output_dir,
                                                                sep=os.path.sep,
                                                                taxid=record['taxonId'])

                fw = None
                if os.path.isfile(fastafile) and not record['taxonId'] in processed_taxon:
                    if action_file_exists == "SA":
                        continue
                    elif action_file_exists == "OA":
                        # this is meant to take care of multiple gbk files in one ncbi dir
                        if d not in processed_dir:
                            fw = open(fastafile, 'w')
                        else:
                            fw = open(fastafile, 'a')
                    else:
                        sys.stdout.write("{file_exists} {taxid}.1.fna".format(file_exists=question_file_exist,
                                                                   taxid=record['taxonId']))
                        action_file_exists = sys.stdin.readline().rstrip()
                        if action_file_exists not in valid_answers_file_exist:
                            action_file_exists = "S"
                        if action_file_exists == 'S' or action_file_exists == 'SA':
                            continue
                        else:
                            fw = open(fastafile, 'w')
                        if action_file_exists == "S" or action_file_exists == "O":
                            action_file_exists = "A"
                else:
                    fw = open(fastafile, 'a')

                if fw is None:
                    sys.stdout.write("Something went wrong in file opening {file}".format(file=fastafile))
                    continue

                # write to file
                l = ">gi{gi}|gb|{acc}|{definition}\n{seq}\n".format(gi=record['gi'],
                                                                  acc=record['accession'],
                                                                  definition=record['definition'],
                                                                  seq=record['seq'])
                fw.write(l)
                fw.close()
                processed_taxon.add(record['taxonID'])
                processed_dir.add(d)

        sys.stdout.write("done\nProcessing of raw NCBI data finished.\n")
