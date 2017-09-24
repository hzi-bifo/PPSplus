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
prediction with already constructed models
arguments: fasta file to predict labels for and configuration file (same one used for training) [model directory]
"""

import os
import sys
import utils
import copy
import argparse
import taxonomy_ncbi
import Bio
import config as c
import parallel


class Predict:
    def __init__(self, cfg):
        """
        :param c: object of Config
        """

        self.workingdir = os.path.dirname(os.path.realpath(sys.argv[0]))
        self.config = cfg
        self.modeldir = ""
        self.models = {}    # mapping fragment length -> model file

        self.sqlite_taxonomy = None

        self.fastafile = None
        self.fastafile_filtered = None

        self.classifier = {}

        self.combined_output_file = ""

    def main_process(self, fastafile, modeldir=None):
        assert os.path.isfile(fastafile) and os.path.isfile(self.config.ncbi_tax_db)

        self.sqlite_taxonomy = taxonomy_ncbi.TaxonomyNcbi(self.config.ncbi_tax_db)

        if modeldir is None:
            modeldir = os.path.join(self.config.settings["project_dir"], "models")
        self.modeldir = modeldir
        self.models = utils.models_from_dir(modeldir, self.config.settings["fragment_len"])

        if self.config.settings["n_classifiers"] > len(self.config.settings["fragment_len"]):
            sys.stdout.write("\nMore classifiers requested than available...changing to {}\n".format(str(len(self.config.settings["fragment_len"]))))
            self.config.settings["n_classifiers"] = len(self.config.settings["fragment_len"])

        self.fastafile = fastafile
        fasta_file_new = "{}.filtered.fna".format(fastafile)
        utils.filter_fastafile(max(self.config.settings["kmer"]), fastafile, fasta_file_new)
        self.fastafile_filtered = fasta_file_new

        self.generate_kmers()
        self.classifier = self.get_classifier()
        self.predict()

        if self.config.settings["postprocessing"]:
            self.post_processing()

        self.sqlite_taxonomy.close()

    def generate_kmers(self):
        kmer_strings = map(lambda x: str(x), self.config.settings["kmer"])
        sys.stdout.write("\nGenerating k-mer features ({})...\n".format("-".join(kmer_strings)))

        fasta2kmers_command = "{script} -s 1 -o 1 -h 1 -l 0 -C {rm} -t {n} -r {rev}".format(script=utils.path_to_binary(self.workingdir, "fasta2kmers"),
                                                                                            rm=self.config.settings["rm_rev_complement"],
                                                                                            n=self.config.settings["kmer_normalization"],
                                                                                            rev=self.config.settings["rev_complement"])
        fasta2kmers2_command = "{script} -a w -s 1 -l 2 -o 1 -b 1 " \
                               "-R {rm} -h 1 " \
                               "-n {n} " \
                               "-r {rev}".format(script=utils.path_to_binary(self.workingdir, "fasta2kmers2"),
                                                  rm=self.config.settings["rm_rev_complement"],
                                                  n=self.config.settings["kmer_normalization"],
                                                  rev=self.config.settings["rev_complement"])

        if len(self.config.settings["kmer"]) == 1:
            k = self.config.settings["kmer"][0]
            fasta2kmers_command = "{cmd} -k {kmer}".format(cmd=fasta2kmers_command, kmer=k)
            fasta2kmers2_command = "{cmd} -k {kmer} -j {kmer}".format(cmd=fasta2kmers2_command, kmer=k)
        else:
            k = max(self.config.settings["kmer"])
            j = min(self.config.settings["kmer"])
            fasta2kmers_command = "{cmd} -k {k} -j {j}".format(cmd=fasta2kmers_command, k=k, j=j)
            fasta2kmers2_command = "{cmd} -k {k} -j {j}".format(cmd=fasta2kmers2_command, k=k, j=j)

        fasta2kmers_command_final = "{cmd} {input} | sed 's/^/1 /' > {output}.sl".format(cmd=fasta2kmers_command,
                                                                                         input=self.fastafile_filtered,
                                                                                         output=self.fastafile)
        fasta2kmers2_command_final = "{cmd} -i {input} -f {output}.sl".format(cmd=fasta2kmers2_command,
                                                                              input=self.fastafile_filtered,
                                                                              output=self.fastafile)

        use_fasta2kmers2 = int(self.config.settings["kmer_normalization"]) <= 1
        if use_fasta2kmers2:
            cmd = fasta2kmers2_command_final
            s = os.system(fasta2kmers2_command_final)
        else:
            cmd = fasta2kmers2_command_final
            s = os.system(fasta2kmers_command_final)
        if s != 0:
            sys.stderr.write("Generating kmers failed with command:\n{cmd}".format(cmd=cmd))
            sys.exit(1)
        sys.stdout.write("done\n")

    def get_classifier(self):
        """
            for each sequence in the fasta file find the closest classifier
                    = the one that was trained for the closest fragment length
            returns a dictionary including index of the classifier (equal to index of the fragment length) as keys
            and a list of indices of the sequences in the fasta file
            that are supposed to be predicted with this classifier.
        """
        fr_kmer = open("{}.sl".format(self.fastafile), "r")
        kmer_lines = fr_kmer.readlines()
        if "" in kmer_lines:    # remove all empty lines
            kmer_lines.remove("")
        fr_kmer.close()

        # go through the fasta file and separate the entries to be classified with different classifiers
        seq_classifier = {}
        n_seq = 0
        for record in Bio.SeqIO.parse(self.fastafile_filtered, 'fasta'):
            seq_len = len(str(record.seq))
            diff_classifier = []
            # get the closet classifier
            for fl in self.config.settings["fragment_len"]:
                diff_classifier.append(abs(fl-seq_len))
            m = min(diff_classifier)
            idx = diff_classifier.index(m)
            # now write to the classifiers file
            fl = self.config.settings["fragment_len"][idx]
            fa = open("{fasta}.{fraglen}.sl".format(fasta=self.fastafile, fraglen=fl), "a")
            # append the corresponding kmer line
            fa.write(kmer_lines[n_seq])
            fa.close()
            n_seq += 1

            if idx in seq_classifier:
                seq_classifier[idx].append(n_seq)
            else:
                seq_classifier[idx] = [n_seq]
        return seq_classifier

    def predict(self):
        sys.stdout.write("Predicting...\n")
        outfiles = set()

        # if commands should be run in parallel, store them in a list for each output file
        # run each list in parallel afterwards
        command_by_outputfile = {}

        for index in self.classifier.keys():
            if len(self.classifier[index]) == 0:
                continue
            self.classify(index, outfiles, command_by_outputfile)

        if self.config.settings["processors"] != 1:
            for commandlist in command_by_outputfile.values():
                if parallel.reportFailedCmd(parallel.runCmdParallel(commandlist, maxProc=self.config.settings["processors"])) is not None:  # Ivan change
                    sys.exit(-1)  # Ivan change
        # join all outputs
        self.combined_output_file = "{}.out".format(self.fastafile)
        utils.concat_files(outfiles, self.combined_output_file)

    def classify(self, index, outfiles, command_by_outputfile):
        classifier_command = "{} -v 0 ".format(utils.path_to_binary(self.workingdir, "svm_phylo_classify"))
        ensemble_command = "{} -v 0 -e 1 ".format(utils.path_to_binary(self.workingdir, "svm_phylo_classify_ensemble"))
        fl = self.config.settings["fragment_len"][index]
        test_file = "{fasta}.{fl}.sl".format(fasta=self.fastafile, fl=fl)
        out_file = "{fasta}.{fl}.out".format(fasta=self.fastafile, fl=fl)
        outfiles.add(out_file)
        # get what classifiers to use
        classifiers_to_use = self.config.settings["fragment_len"][index:min(index + self.config.settings["n_classifiers"],
                                                                            len(self.config.settings["fragment_len"]))]

        classifier_to_use_string = ""
        for c in classifiers_to_use:
            classifier_to_use_string += str(c)
            classifier_to_use_string += ","
        classifier_to_use_string = classifier_to_use_string[:len(classifier_to_use_string)-1]
        sys.stdout.write("\tfragments close to length {fl} with classifiers {c}\n".format(fl=fl,
                                                                                          c=classifier_to_use_string))
        n_classifiers_to_use = len(classifiers_to_use)
        if n_classifiers_to_use == 1:
            extra_command = "{test} {model} {out}".format(test=test_file, model=self.models[fl], out=out_file)
            final_command = "{classifier}{extra}".format(classifier=classifier_command, extra=extra_command)
        else:
            extra_command = "-m {n} {test} ".format(n=n_classifiers_to_use, test=test_file)
            for i in classifiers_to_use:
                extra_command = "{cmd}{model} ".format(cmd=extra_command, model=self.models[i])
            extra_command = "{cmd}{out} ".format(cmd=extra_command, out=out_file)
            final_command = "{cmd}{ext}".format(cmd=ensemble_command, ext=extra_command)

        if self.config.settings["processors"] != 1:
            sys.stdout.write("running predictions in parallel\n")
            if out_file not in command_by_outputfile:
                command_by_outputfile[out_file] = [parallel.TaskCmd(final_command)]
            else:
                command_by_outputfile[out_file].append(parallel.TaskCmd(final_command))
        else:
            s = os.system(final_command)
            if s != 0:
                sys.stderr.write("Error in classification: couldn't run the system command.\n")
                sys.exit(1)

    def post_processing(self):
        sys.stdout.write("Post-processing... \n")
        fr = open(self.combined_output_file, "r")
        paths = {}
        preds = []
        heads = []
        for line in fr:
            if line == "":
                continue
            line = line.rstrip().split("\t")
            heads.append(line[0])
            preds.append(line[4])
            if line[4] not in paths:
                for rank in taxonomy_ncbi.TAXONOMIC_RANKS:
                    parent_at_rank = self.sqlite_taxonomy.parentAtRank(line[4], rank)
                    if parent_at_rank is None:
                        parent_at_rank = ""
                    paths[preds].append(parent_at_rank)
        fr.close()
        # create a PP style output file
        out_file_pp = "{}.PP.out".format(self.fastafile)
        # write the file
        # after converting taxonomy ids to names
        fw = open(out_file_pp, "w")
        taxonomy_ncbi.TAXONOMIC_RANKS.reverse()
        fw.write("#Output in PhyloPythia format\n"
                 "#Taxonomic Hierarchy: NCBI Taxonomy\n"
                 "#Prediction method: PhyloPythiaS\n"
                 "#Contact: PhyloPythiaS mailing list (phylopythias@uni-duesseldorf.de)\n"
                 "#ID\t{}".format("\t".join(taxonomy_ncbi.TAXONOMIC_RANKS)))
        taxonomy_ncbi.TAXONOMIC_RANKS.reverse()

        paths_name = {}
        for i in range(len(preds)):
            pred = preds[i]
            head = heads[i]
            if pred not in paths_name:
                path = paths[pred]
                path_name = []
                for p in path:
                    n = self.sqlite_taxonomy.getScientificName(p)
                    if n is None or n == "":
                        n = "\t"
                    path_name.append(n)
                path_name.reverse()
                paths_name[pred] = path_name
            else:
                path_name = paths_name[pred]
            path_for_output = copy.copy(path_name)
            path_for_output.append(head)
            path_for_output.reverse()
            fw.write("\t".join(path_for_output))
        fw.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="configuration file same as used for training", action='store', required=True)
    parser.add_argument("-fasta", help="fastafile including sequences, that should be predicted",
                        action='store', required=True)
    parser.add_argument("-models", help="directory including models", action='store', default=None)
    args = parser.parse_args()

    p = Predict(c.Config(args.c))
    p.main_process(args.fasta, args.models)
