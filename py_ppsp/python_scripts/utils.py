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

About: utilities for handling files, folders and data structures

"""
import sys
import os
import time
import re
import Bio
from Bio import SeqIO
from types import *


def get_answer_timeout(question=""):
    sys.stdout.write(question)
    start_time = time.time()  # this is time in seconds
    while True:
        answer = sys.stdin.readline().strip().lower()
        if answer == "y":
            return "y"
        curr_time = time.time()
        time_ellapsed = curr_time - start_time
        if time_ellapsed is 120 or answer == "n":
            return "n"


def any_list_to_string(l, sep="\t"):
    # unit tests
    try:
        s = sep.join(l)
    except TypeError:
        s = ""
        for e in l:
            if str(e) == "":
                continue
            s += str(e)
            s += sep
        s = s[:-1]
    return s


def get_lines(f):
    """
        get all lines from a file

        (unit tests created)
    """
    try:
        fr = open(f)
    except:
        sys.stderr.write("can not open file {}\n".format(f))
        sys.exit(1)
    lines = fr.readlines()
    fr.close()
    return lines


def dir_is_empty(d):
    """
        check if a directory is empty

        (unit tests created)
    """
    assert os.path.isdir(d)

    entries = os.listdir(d)
    if len(entries) is 0:
        return True
    else:
        return False


def list_to_file(l, filename):
    # unit tests included
    try:
        fw = open(filename, "w")
    except:
        sys.stdout.write("error writing file {}\n".format(filename))
        return False

    s = any_list_to_string(l, os.linesep)
    fw.write(s)

    fw.close()
    return True


def filter_sequence(kmer_max, seq):
    """
        replaces non-ACGT in a given sequence by XXX
    :param kmer_max: number of X's
    :param seq: sequence that should be filtered
    :return: filtered sequence
    """
    assert type(kmer_max) is IntType, "kmer_max must be integer"
    assert (kmer_max > 0), "kmer_max must be > 0"
    assert type(seq) is StringType, "DNA sequences should be strings"
    assert len(seq) > 0, "the sequence should not be empty"

    replace = "X" * kmer_max
    pattern = "[^ACGT]{" + str(kmer_max) + ",}"
    pat = re.compile(pattern, flags=re.IGNORECASE)
    return re.sub(pat, replace, seq)


def filter_fastafile(kmer_max, fasta_original, fasta_new):
    """
        for all sequences in a given fastafile:
            replaces non-ACGT in sequence by XXX
    :param kmer_max: number of X's
    """
    assert type(kmer_max) is IntType, "kmer_max must be integer"
    assert (kmer_max > 0), "kmer_max must be > 0"

    handle = open(fasta_new, "w")
    records = []
    for sequence in Bio.SeqIO.parse(fasta_original, 'fasta'):
        sequence_string = sequence.seq.tostring()
        sequence_new = filter_sequence(kmer_max, sequence_string)
        seq_id = sequence.id
        records.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(sequence_new), id=seq_id, description=""))
    Bio.SeqIO.write(records, handle, "fasta")
    handle.close()


def concat_files(filenames, outputfile):
    with open(outputfile, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def fasta2kmers_base(config, workingdir):
    fasta2kmers_command = "{bin} " \
                          "-s 1 -o 1 -h 1 -l 0 " \
                          "-C {rm_rev_comp} " \
                          "-t {norm} " \
                          "-r {rev_comp}".format(bin=path_to_binary(workingdir, "fasta2kmers"),
                                                 rm_rev_comp=config.settings["rm_rev_complement"],
                                                 norm=config.settings["kmer_normalization"],
                                                 rev_comp=config.settings["rev_complement"])
    # identify what kmers
    if len(config.settings["kmer"]) == 1:
        k = config.settings["kmer"][0]
        fasta2kmers_command = "{cmd} -k {kmer}".format(cmd=fasta2kmers_command, kmer=k)
    else:
        k = max(config.settings["kmer"])
        j = min(config.settings["kmer"])
        fasta2kmers_command = "{cmd} -k {k} -j {j}".format(cmd=fasta2kmers_command, k=k, j=j)
    return fasta2kmers_command


def fasta2kmers2_base(config, workingdir):
    fasta2kmers2_command = "{bin} -a a -s 1 -l 1 -o 1 -b 1 " \
                           "-R {rm_rev_comp} " \
                           "-n {norm} " \
                           "-r {rev_comp}".format(bin=path_to_binary(workingdir, "fasta2kmers2"),
                                                  rm_rev_comp=config.settings["rm_rev_complement"],
                                                  norm=config.settings["kmer_normalization"],
                                                  rev_comp=config.settings["rev_complement"])
    # identify what kmers
    if len(config.settings["kmer"]) == 1:
        k = config.settings["kmer"][0]
        fasta2kmers2_command = "{cmd} -k {kmer} -j {kmer}".format(cmd=fasta2kmers2_command, kmer=k)
    else:
        k = max(config.settings["kmer"])
        j = min(config.settings["kmer"])
        fasta2kmers2_command = "{cmd} -k {k} -j {j}".format(cmd=fasta2kmers2_command, k=k, j=j)

    return fasta2kmers2_command


def path_to_binary(workingdir, binary_name):
    # unit tests added
    return os.path.join(workingdir, os.path.pardir, "bin", binary_name)


def models_from_dir(modeldir, frag_lens):
    # added unit tests

    assert ListType(frag_lens) and len(frag_lens) > 0
    assert os.path.isdir(modeldir)

    m = {}
    # check if all the model files exist in the project directory
    model_files = os.listdir(modeldir)
    for fl in frag_lens:
        for mf in model_files:
            mf = os.path.basename(mf)

            if "_" not in mf:
                continue

            fl_mf = int(mf.split("_")[0])
            if int(fl_mf) == int(fl):
                m[fl] = os.path.join(modeldir, mf)
        try:
            m[fl]
        except IndexError:
            sys.stderr.write("could not find model for fragment length: {}\n".format(fl))
            sys.exit(1)
    return m


def window(seq_concat, winsize, indices, stepsize=1):
    """
            get fragments from a given concatenated sequence at given indices
            by specifying a shorter step size than window size you will get overlapping fragments
            e.g. if the given sequence can be divided into
                frag0 frag1 frag2 frag3 frag4
                and the given indices are 1 and 4, it will return frag1 and frag4
            return:
                (start position of the fragment in the sequence, fragment)

            (added unit tests)
    """
    assert winsize > 0
    assert len(indices) > 0
    assert stepsize > 0

    results = []
    for index in indices:
        start = index * stepsize
        stop = start + winsize
        if stop > len(seq_concat):
            break
        window = seq_concat[start : stop]
        results.append((start, window))
    return results


def get_sequence_infos(fastafile):
    """
        input:
            fastafile with one or more sequences:
                >sequence name 1
                    sequence 1
                 ...
                >sequence name n
                    sequence n
        return:
            concatenated sequence: sequence1 N sequence2 N ... sequence n
            a name for the concatenated sequence (= definition)
    """
    assert os.path.isfile(fastafile)

    sequences = []
    definition = ""
    for sequence in Bio.SeqIO.parse(fastafile, 'fasta'):
        sequences.append(sequence.seq.tostring())
        if definition == "":
            definition = sequence.id
    seq_concat = "N".join(sequences)
    seq_concat = re.sub('[^a-zA-Z]', '', seq_concat)
    return seq_concat, definition


def write_fragments(fl, step, seq_concat, definition, mapped, fastafile, sampled_frag):
    
    if not (len(sampled_frag) > 0):  # Ivan added
        # print 'write_fragments:', fl, step, len(seq_concat), definition, mapped, fastafile, sampled_frag  # Ivan added
        # print('Utils, write_fragments: problem sampled_frag: %s' % len(sampled_frag))  # Ivan added
        if fl <= len(seq_concat):
            sampled_frag = [0]  # Ivan added
            # print('write_fragments: ok')
        else:
            print('write_fragments: skipped')
            print 'write_fragments:', fl, step, len(seq_concat), definition, mapped, fastafile, sampled_frag  # Ivan added
            print('Utils, write_fragments: problem sampled_frag: %s' % len(sampled_frag))  # Ivan added
            return  # Ivan added

    assert IntType(fl) and fl > 0
    assert IntType(step) and step > 0
    assert ListType(sampled_frag) and len(sampled_frag) > 0
    assert StringType(seq_concat) and len(seq_concat) > 0

    if definition is None or definition == "<unknown description>":
        definition = ""
    if mapped is None:
        mapped = ""

    records = []
    for start, s in window(seq_concat, fl, sampled_frag, stepsize=step):
        if s == "":
            continue
        rec = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s), id=definition+"|at"+str(start)+"|label:"+mapped, description="")
        lower_rec = rec.lower()
        records.append(lower_rec)

    if os.path.isfile(fastafile):
        handle = open(fastafile, "a")
    else:
        handle = open(fastafile, "w")

    writer = Bio.SeqIO.FastaIO.FastaWriter(handle, wrap=1000)
    writer.write_file(records)
    handle.close()


"""
        Unit tests for the different functions
"""

import unittest
import random
import tempfile
import shutil


class Testing(unittest.TestCase):

    def test_any_list_to_string(self):
        listings = [[1, 2, 3, 4, 5, 6],
                    ["1", "2", "3", "4", "5", "6"],
                    ["1", "2", 3, "4", 5, "6"],
                    ["1", "2", "3", "4", "", "5", "6"],
                    [1, 2, 3, 4, "", 5, 6]]
        s = "1,2,3,4,5,6"
        for l in listings:
            self.assertEqual(any_list_to_string(l, ","), s)

    def test_get_lines(self):
        l = []
        for i in range(100):
            random_number = random.randint(0, 10000000)
            l.append(str(random_number)+"\n")
        fw = open("test.txt", "w")
        fw.write(l)
        fw.close()
        read_lines = get_lines(fw)
        self.assertListEqual(read_lines, l)
        os.remove("test.txt")

    def test_dir_is_empty(self):
        directory = tempfile.mkdtemp()
        # this directory should be empty
        self.assertEqual(True, dir_is_empty(directory))
        # create an empty file
        fw = open(os.path.join(directory, "test.txt"), "w")
        fw.close()
        self.assertEqual(False, dir_is_empty(directory))
        os.remove(os.path.join(directory, "test.txt"))
        # create an empty directory
        os.mkdir(os.path.join(directory, "hanny"))
        self.assertEqual(False, dir_is_empty(directory))
        shutil.rmtree(directory)

    def test_list_to_file(self):
        listings = [[1, 2, 3, 4, 5, 6],
                    ["1", "2", "3", "4", "5", "6"],
                    ["1", "2", 3, "4", 5, "6"],
                    ["1", "2", "3", "4", "", "5", "6"],
                    [1, 2, 3, 4, "", 5, 6]]
        lines = ["1" + os.linesep,
                 "2" + os.linesep,
                 "3" + os.linesep,
                 "4" + os.linesep,
                 "5" + os.linesep,
                 "6" + os.linesep]
        for l in listings:
            list_to_file(l, "test.txt")
            with open("test.txt") as input:
                read_lines = input.readlines()
            self.assertListEqual(lines, read_lines)
            os.remove("test.txt")

    def test_path_to_binary(self):
        directory = tempfile.mkdtemp()

        workingdir = os.path.join(directory, "scripts")
        os.mkdir(workingdir)

        bin_dir = os.path.join(directory, "bin")
        os.mkdir(bin_dir)

        text = "this is a binary file, that has been called."
        with open(os.path.join(bin_dir, "test.bin")) as outfile:
            outfile.write(text)

        with open(path_to_binary(workingdir, "test.bin")) as inputfile:
            l = inputfile.readline().strip()
            self.assertEqual(text, l)

        shutil.rmtree(directory)

    def test_models_from_dir(self):
        directory = tempfile.mkdtemp()
        fl_to_file = {}
        frag_lens = [100, 200, 300, 400]
        for fl in frag_lens:
            modelfile = os.path.join(directory), "{}_10.sl".format(str(fl))
            fl_to_file[fl] = modelfile
            fw = open(modelfile, "w")
            fw.close()
        self.assertDictEqual(fl_to_file, models_from_dir(directory, frag_lens))

        # missing model for one fragment len -> exit
        frag_lens.append(10000)
        with self.assertRaises(SystemExit) as cm:
            models_from_dir(directory, frag_lens)
        self.assertEqual(cm.exception.code, 1)

        with self.assertRaises(AssertionError):
            models_from_dir(directory, [])
            models_from_dir(directory, False)
            models_from_dir(directory, 1)
            models_from_dir(directory, "hallo?")
            models_from_dir("i am not a directory", [100, 200, 300])

    def test_window_non_overlapping(self):
        seq = "ACGTTTTGTACTGTACGTTTTTGACGTCCCCCCGT"
        # non overlapping fragments
        # all possible fragments results[startpos]=fragment
        results = {0: "ACGTT",
                   5: "TTGTA",
                   10: "CTGTA",
                   15: "CGTTT",
                   20: "TTGAC",
                   25: "GTCCC",
                   30: "CCCGT"}
        ct = 0
        for startindex, frag in window(seq, 5, range(7), 5):
            ct += 0
            self.assertEqual(frag, results[startindex])
        self.assertEqual(ct, 7)  # found all 7 fragments?

        # only fragments 2, 3 and 6 (because we need only 3)
        results = {10: "CTGTA",
                   15: "CGTTT",
                   30: "CCCGT"}
        ct = 0
        for startindex, frag in window(seq, 5, [2, 3, 6], 5):
            ct += 1
            self.assertEqual(frag, results[startindex])
        self.assertEqual(ct, 3)     # found all 3 fragments?

    def test_window_overlapping_step1(self):
        seq = "ACGTTTTGTACTGTACGT"
        # overlapping fragments
        # all possible fragments results[startpos]=fragment
        results = {0: "ACGTTTTGTA",
                   1: "CGTTTTGTAC",
                   2: "GTTTTGTACT",
                   3: "TTTTGTACTG",
                   4: "TTTGTACTGT",
                   5: "TTGTACTGTA",
                   6: "TGTACTGTAC",
                   7: "GTACTGTACG",
                   8: "TACTGTACGT"}
        ct = 0
        for startindex, frag in window(seq, 5, range(10), 1):
            ct += 0
            self.assertEqual(frag, results[startindex])
        self.assertEqual(ct, 10)  # found all 10 fragments?

    if __name__ == '__main__':
        unittest.main()
