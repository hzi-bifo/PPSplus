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
generate training data and train models on it
input file is the configuration file
"""

import sys
import os
import re
import random
import glob
import argparse
import utils
import taxonomy_ncbi
import config
import parallel
import math
from status import Status
import process_ncbi
import shutil
import ncbi2newick
from types import *
import logging
import logmethods


class Train:

    def __init__(self, cfg, default_ans=False, backup=False, logged=False):
        self.workingdir = os.path.dirname(os.path.realpath(sys.argv[0]))

        # some mappings
        self.tree_organism_map = {}     # node --> list of organisms which were mapped to that node
        self.organism_tree_map = {}     # organism --> closest node on the tree
        self.organism_file_map = {}     # organism --> list of files including sequences

        self.organisms = set()      # all organisms available
        self.nodes = []             # selected organisms = nodes on the taxonomic tree

        self.genomes_excluded = []  # genomes that have not been considered during computation
        self.organisms_invalid = set()  # organisms with a lack of mapping // or that have been mapped to the root

        self.sqlite_taxonomy = None     # SQLite database including the NCBI taxonomy
        self.config = cfg  # contains all settings

        # # # # # # # Advanced training options # # # # # # # # # # # # # #
        # learning
        self.loss_function = 1  # 0:0/1 loss, 1:path loss
        self.z_standardization = 1  # 0:no standardization, 1:z standardization
        self.misc_nodes = 1     # 0:no misc nodes, 1:add misc nodes
        self.n_frags_per_node = 0.0     # number of fragments per node

        self.yes = default_ans    # say yes to every question?

        self.tree_file = ""

        self.backup = backup
        self.stat = None
        self.backupdir = None

        if not logged:
            # this will write a file to cwd
            logmethods.initialise('train.log')
        self.log = logging.getLogger('train')   # better root ?

    def from_backup(self, backupdir):
        assert os.path.isdir(backupdir)
        my_log = logging.getLogger('backup')

        if self.backup:
            my_log.warning('The backup will be produced but can not be used for another run.')

        self.stat = Status(logged=True)
        try:
            self.stat.load_backup(backupdir)
        except SystemExit:
            my_log.critical('Could not load backup from:\n{}'.format(backupdir))
            sys.exit(1)

        pipeline_position = self.stat.variables["status"]

        if pipeline_position == 6:
            my_log.critical("Please re-run the pipeline.. "
                             "There is no secure way to find the correct point to start from")
            sys.exit(1)

        # init all variables
        for variable_name in self.stat.variables:
            val = self.stat.variables[variable_name]
            setattr(self, variable_name, val)   # same as self.variable_name = val

        self.sqlite_taxonomy = taxonomy_ncbi.TaxonomyNcbi(self.config.ncbi_tax_db)

        if pipeline_position < 1:
            self.process_ncbi()
        if pipeline_position < 2:
            self.tree_process()
        if pipeline_position < 3:
            self.map_genomes_on_tree()
        if pipeline_position < 4:
            self.generate_seq(done_files=self.stat.files.keys())
        if pipeline_position < 5:
            self.sample_specific(done_files=self.stat.files.keys())
        if pipeline_position < 6:
            self.generate_kmer_features()
        if pipeline_position < 7:
            self.build_models()

    def main_processing(self):
        """
            all steps of the training pipeline
                - processing ncbi sequences
                - tree processing (copying newick tree or building tree from a clade list)
                - mapping genomes on this tree
                - generating fragments from the sequences
                - dealing with sample specific data
                - generating kmer features
                - building models
        """

        if not self.config.settings["only_models"]:
            # if you do not want to build only models, check if the project directory is empty
            if len(os.listdir(self.config.settings["project_dir"])) != 0:
                self.log.warning("The project directory is not empty, this can result in unpredictable behavior.")
                if self.yes:
                    answer = "y"
                else:
                    answer = utils.get_answer_timeout("Remove? [Y/N]")
                if answer == "y":
                    self.log.info("Deleting old project directory..")
                    shutil.rmtree(self.config.settings["project_dir"])
                    os.mkdir(self.config.settings["project_dir"])
                else:
                    self.log.critical("Please provide an empty project directory. Quiting...")
                    sys.exit(1)

            if self.backup:
                self.stat = Status(logged=True)
                self.backupdir = os.path.join(self.config.settings["project_dir"], "backup")
                os.mkdir(self.backupdir)

            self.log.info("creating folder structure...")
            self.create_folderstructure()

            self.log.info("checking database")
            self.check_db()

            if len(self.config.genomes_exclude) != 0:
                # write organisms that will not be considered into a file
                utils.list_to_file(self.config.genomes_exclude,
                             os.path.join(self.config.settings["project_dir"], "excluded.txt"))

                if self.stat is not None:
                    self.stat.add_written_file(os.path.join(self.config.settings["project_dir"], "excluded.txt"))
                    self.stat.write_backup(self.backupdir)

            self.log.info("Processing NCBI data...")
            self.process_ncbi()
            if self.stat is not None:
                self.stat.change_variable(1, "status")
                self.stat.change_variable(self.organism_file_map, "organism_file_map")
                self.stat.change_variable(self.genomes_excluded, "genomes_excluded")
                self.stat.change_variable(self.organisms, "organisms")
                self.stat.write_backup(self.backupdir)

            self.log.info("tree processing...")
            self.tree_process()
            if self.stat is not None:
                self.stat.change_variable(2, "status")
                self.stat.change_variable(self.nodes, "nodes")
                self.stat.change_variable(self.tree_file, "tree_file")
                self.stat.write_backup(self.backupdir)

            self.log.info("mapping genomes on the tree...")
            self.map_genomes_on_tree()
            if self.stat is not None:
                self.stat.change_variable(3, "status")
                self.stat.change_variable("n_frags_per_node", "status")
                self.stat.change_variable("tree_organism_map", "status")
                self.stat.change_variable("organisms_invalid", "status")
                self.stat.change_variable("organism_tree_map", "status")
                self.stat.write_backup(self.backupdir)

            self.log.info("generating sequence fragments...")
            self.generate_seq()
            if self.stat is not None:
                self.stat.change_variable(4, "status")
                self.stat.write_backup(self.backupdir)

            self.log.info("sample specific stuff...")
            self.sample_specific()
            if self.stat is not None:
                self.stat.change_variable(5, "status")
                self.stat.write_backup(self.backupdir)

            self.log.info("generating kmer features...")
            self.generate_kmer_features_concat()
            if self.stat is not None:
                self.stat.change_variable(6, "status")
                self.stat.write_backup(self.backupdir)

        else:
            self.log.info("reading tree string")
            self.config.settings["tree_file"] = os.path.join(self.config.settings["project_dir"],
                                                             "tree.newick")

        self.log.info("building models")
        self.build_models()

        if self.config.settings["clean_up_train"]:
            self.log.info("Cleaning..")
            shutil.rmtree(os.path.join(self.config.settings["project_dir"], "train_data"))
            shutil.rmtree(os.path.join(self.config.settings["project_dir"], "sampled_fasta"))

        self.sqlite_taxonomy.close()
        self.log.info("Processing finished ...models are ready in {}".format(os.path.join(self.config.settings["project_dir"], "models")))

    def process_ncbi(self):
        """
            download NCBI sequences,
            parse gbk files to fasta format,
            rename sequences and sort per genome
        """
        my_log = logging.getLogger('train:process_ncbi')
        if utils.dir_is_empty(self.config.settings["ncbi_processed_dir"]):
            my_log.warning("The NCBI_PROCESSED_DIR is empty, it is possible to download data from NCBI.")
            my_log.info("I can automatically download Bacterial & Archael data "
                        "(for more possibilities see INSTALL.txt).")

            if self.yes:
                answ = "y"
            else:
                my_log.info("Download sequence data from NCBI? [Y/N] (default=Y, timeout 2 minutes)")
                answ = utils.get_answer_timeout()

            if answ != "y":
                error_message = "There is no training data available, provide some and run the program again." \
                                "Read INSTALL.txt for details on how this can be done."
                my_log.critical(error_message)
                sys.exit(1)

            my_log.info("Download may take some time ...")
            os.mkdir(os.path.join(self.config.settings["project_dir"], "tmp"))

            tmp_dir = os.path.join(self.config.settings["project_dir"], 'tmp')
            data_archive = os.path.join(tmp_dir, 'all.gbk.tar.gz')

            success = os.system("wget -O {} ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gbk.tar.gz".format(data_archive))
            if success != 0:
                my_log.critical("Error in downloading sequence data from NCBI.")
                sys.exit(1)
                # unpack data

            unpack_cmd = "tar xfz {a} -C {tmp}".format(a=data_archive, tmp=tmp_dir)
            success = os.system(unpack_cmd)

            if success != 0:
                my_log.critical("Error in unpacking the downloaded sequence data.")
                sys.exit(1)

            # process the data and create the fasta files in ncbi_dir
            process_object = process_ncbi.Process_NCBI(tmp_dir, self.config.settings["ncbi_processed_dir"])
            success = process_object.run()
            if not success:
                sys.exit(1)

            # clean the dowloaded NCBI data
            shutil.rmtree(tmp_dir)

        # get all the organism names from the files in the ncbi_dir
        # this can be used for generating generic clades
        n_sequences = 0
        files = glob.glob("{dir}{sep}*.*".format(dir=self.config.settings["ncbi_processed_dir"],
                                                 sep=os.path.sep))
        for f in files:
            ext = f.split(".")[-1]
            if len(self.config.settings["extensions"]) > 0 or \
                    (len(self.config.settings["extensions"]) == 1 and len(self.config.settings["extensions"][0]) > 0):
                        # extensions entweder [] oder [""]
                if ext is not None and ext not in self.config.settings["extensions"]:
                    continue
            if "." not in f:
                my_log.debug("Invalid file: {}..skipping".format(f))
                continue
            else:
                organism = f.split(os.path.sep)[-1].split(".")[0]

            # exclude this genome if asked to
            if organism in self.config.genomes_exclude:
                self.genomes_excluded.append(organism)
                continue

            n_sequences += 1
            self.organisms.add(organism)
            if organism not in self.organism_file_map:
                self.organism_file_map[organism] = [f]
            else:
                self.organism_file_map.append(f)

        if len(self.genomes_excluded) is not 0:
            my_log.info("excluded {} sequences.".format(str(len(self.genomes_excluded))))

    def tree_process(self):
        my_log = logging.getLogger('train:tree_processing')
        if self.config.settings["tree_file"] is None or self.config.settings["tree_file"] == "":
            # create clade list (own method?)
            descandents = {}
            for organism in self.organisms:
                for rank in self.config.settings["taxonomy_ranks"]:
                    parent = self.sqlite_taxonomy.parent_at_rank(organism, rank)
                    print('!!! Call: parent = self.sqlite_taxonomy.parent_at_rank(organism, rank) !!!')
                    if parent in descandents:
                        descandents[parent] += 1
                    else:
                        descandents[parent] = 1
            for parent in descandents.keys():
                if parent is None or parent == "":
                    continue
                if descandents[parent] >= self.config.settings["n_min_genomes_generic"]:
                    self.nodes.append(parent)
        else:
            # read the tree_string
            nodes_tmp = utils.get_lines(self.config.settings["tree_file"])
            for n in nodes_tmp:
                if n != "" and n is not None:
                    self.nodes.append(n.rstrip("\n"))
            tree_string = self.nodes[0]

        # check if the tree is a newick string or a node list
        t_file = os.path.join(self.config.settings["project_dir"], "tree.newick")
        if ";" not in self.nodes[0] and len(self.nodes) >= 1:
            my_log.debug("Generating tree from the clades list ({} clades)....".format(str(len(self.nodes))))
            clades_file = os.path.join(self.config.settings["project_dir"], "clades.txt")
            utils.list_to_file(self.nodes, clades_file)

            if self.stat is not None:
                self.stat.add_written_file(clades_file)
                self.stat.succesfully_written(clades_file)

            # run script to create tree
            obj = ncbi2newick.Ncbi2Newick(self.config.ncbi_tax_db, logged=True)
            obj.tree_from_nodes(clades_file)
            obj.tree_to_file(t_file)
            obj.close()

        else:
            my_log.debug("Copying tree to the project directory...")
            fw = open(t_file, "w")
            if self.stat is not None:
                self.stat.add_written_file(t_file)
            fw.write(tree_string)
            if self.stat is not None:
                self.succesfully_written(t_file)
            fw.close()

            # change tree_file to the tree in the project directory
        # get back the tree_string  and nodes, this is necessary for further processing
        self.tree_file = os.path.join(self.config.settings["project_dir"], "tree.newick")
        fr = open(self.tree_file, "r")
        tree_string = fr.readline().rstrip()
        fr.close()
        if tree_string == "":
            my_log.critical("First line is empty in the newick file: {}".format(self.tree_file))
            sys.exit(1)
        self.nodes = ncbi2newick.get_nodes_from_newick(self.tree_file)

    def map_genomes_on_tree(self):
        my_log = logging.getLogger('train:map_genomes')

        self.organism_tree_map = {}
        mappings = set()

        for organism in self.organisms:
            if organism in self.organism_tree_map:
                continue    # already exists
            mapped = self.get_mapped_organism(organism)
            if mapped is None:
                self.organisms_invalid.add(organism)
                continue
            if str(mapped) == "1":
                self.organisms_invalid.add(organism)

            self.organism_tree_map[organism] = mapped
            if mapped not in self.tree_organism_map:
                self.tree_organism_map[mapped] = set()
            self.tree_organism_map[mapped].add(organism)
            mappings.add(mapped)

        self.n_frags_per_node = int(math.ceil(self.config.settings["number_examples"]/len(mappings)))

        if len(self.organisms_invalid) >= 1:
            my_log.info("these {nr} organisms will not be processed due to lack of mapping:\n{orgs}"
                             .format(nr=str(len(self.organisms_invalid)), orgs="\t".join(self.organisms_invalid)))
            utils.list_to_file(self.organisms_invalid,
                               os.path.join(self.config.settings["project_dir"], "organisms_not_used.txt"))
            if self.stat is not None:
                self.stat.add_written_file(os.path.join(self.config.settings["project_dir"], "organisms_not_used.txt"))
                self.stat.succesfully_written(os.path.join(self.config.settings["project_dir"], "organisms_not_used.txt"))

    def generate_seq(self, done_files=[]):
        my_log = logging.getLogger('train:generate_seq')
        # tree_organism_map contains node as keys and organism array as values
        nodesdone = 0
        for node in self.tree_organism_map:
            orgs = self.tree_organism_map[node]
            files = []
            pairs = []
            for o in orgs:
                if o in self.organisms_invalid:
                    # invalid organism
                    continue
                for f in self.organism_file_map[o]:
                    files.append(f)
                    pairs.append((os.path.getsize(f), f))

            if len(files) == 0:
                my_log.debug("no files for this node")
                continue

            n_frag_per_file = max(int(math.ceil(self.n_frags_per_node/len(files))), 1)
            if n_frag_per_file == 0:
                continue

            nodesdone += 1
            my_log.info("processing node: {lab} ({current}/{of})".format(current=str(nodesdone),
                                                                         lab=node,
                                                                         of=str(len(self.tree_organism_map.keys()))))
            my_log.debug("files for this node:{nr}\tfragments per file:{n_frag}".format(nr=str(len(files)),
                                                                                        n_frag=str(n_frag_per_file)))

         #   frag_len_to_frags_so_far = {}  
            # how many fragments have been taken for each fragment length?
            if n_frag_per_file == 1:  # Ivan added
                # sort list of files by file size
                # this means you will get fragments for organisms with much data first
                pairs.sort(key=lambda s: s[0], reverse=True)  # Ivan added
                pairs = pairs[:self.n_frags_per_node]  # Ivan added
                files = []  # Ivan added
                for p in pairs:  # Ivan added
                    files.append(p[1])  # Ivan added

            for f in files:
                filename = os.path.basename(f)
                # read in fasta sequences and definitions from ncbi_dir/.....
                seq_concat, definition = utils.get_sequence_infos(f)
                #     remove non-ACGT to get a better normalized k-mer vector
                seq_concat = utils.filter_sequence(max(self.config.settings["kmer"]), seq_concat)

                for fl in self.config.settings["fragment_len"]:
                    if len(seq_concat) < fl:
                        continue
                    # determine how many fragments to take

                    n_frag = int(math.floor(len(seq_concat)/fl))
                    if n_frag == 0:
                        continue

                    sampled_frag = range(n_frag)
                    if n_frag > n_frag_per_file:
                        # choose some random fragments
                        # get randomized list sampled_frag
                        # set the seed to sequence length to make the functionality deterministic
                        random.seed(len(seq_concat))
                        random.shuffle(sampled_frag)
                        # choose n_frag_per_file random elements
                        sampled_frag = sampled_frag[0:n_frag_per_file]

                    sample_dir = os.path.join(self.config.settings["project_dir"], "sampled_fasta")
                    fl_dir = os.path.join(sample_dir, str(fl))

                    fastafile = os.path.join(fl_dir, filename)

                    if fastafile in done_files:
                        continue

              #      if n_frag_per_file == 1:
              #          if fl not in frag_len_to_frags_so_far:
              #              frag_len_to_frags_so_far[fl] = len(sampled_frag)
              #          else:
              #              frag_len_to_frags_so_far[fl] += len(sampled_frag)
              #          if frag_len_to_frags_so_far[fl] > self.n_frags_per_node:  # already took enough fragments
              #              break

                    utils.write_fragments(fl, fl, seq_concat, definition, node, fastafile, sampled_frag)


                # this is very inefficient with the current implementation
                #    if self.stat is not None:
                #        self.stat.add_written_file(fastafile)
                #        self.stat.succesfully_written(fastafile)
                #        self.stat.write_backup(self.backupdir)

    def sample_specific(self, done_files=[]):
        my_log = logging.getLogger('train:sample_specific')
        #     check if sample specific data is provided:
        if self.config.settings["sample_specific_dir"] == "":
            my_log.debug('no sample specific data')
            return
        if utils.dir_is_empty(self.config.settings["sample_specific_dir"]):
            my_log.debug('no sample specific data')
            return

        genomes_exclude_ss = []
        n_sequences_ss = 0
        my_log.info("Processing sample specific data (all data will be used)...")

        for e in os.listdir(self.config.settings["sample_specific_dir"]):
            # valid extension?
            if len(self.config.settings["extensions"]) != 0:
                ext = e.split(".")[-1]
                if ext == "" or ext not in self.config.settings["extensions"]:
                    continue
            organism = e.split(".", 1)[0]
            if organism == "":
                my_log.warning("Invalid sample specific file: {} skipping..".format(e))
                continue

            # map organism
            if organism in self.nodes:
                node = organism
            else:
                node = self.get_mapped_organism(organism)

            if node is None:
                my_log.info("Could not map {} on the tree".format(organism))
                continue
            elif str(node) == "1":
                my_log.info("Skipping {} due to lack of mapping".format(organism))
                continue

            seq_concat, definition = utils.get_sequence_infos(os.path.join(self.config.settings["sample_specific_dir"], e))

            for i in range(len(self.config.settings["fragment_len"])):
                fl = self.config.settings["fragment_len"][i]
                try:
                    step = self.config.settings["sample_specific_step"][i]
                except IndexError:  # shgould not happen at all
                    step = fl

                if step == 0 or step is None:
                    step = fl

                if len(seq_concat) < fl:
                    my_log.debug("No sample specific data for organism {o} at frag_len {fl}".format(o=organism, fl=fl))
                    continue

                number_frags = (len(seq_concat)-fl)/step
                sample_dir = os.path.join(self.config.settings["project_dir"], "sampled_fasta")
                fl_dir = os.path.join(sample_dir, str(fl))
                fastafile = os.path.join(fl_dir, e)

                utils.write_fragments(fl, step, seq_concat, definition, node, fastafile, range(number_frags))

            #    if self.stat is not None:      #  this is very inefficient at the moment !!!
            #        self.stat.succesfully_written(fastafile)
            #        self.stat.write_backup(self.backupdir)

            n_sequences_ss += 1

        if n_sequences_ss == 0:
            my_log.error("no data processed in SAMPLE_SPECIFIC_DIR")
        if len(genomes_exclude_ss) != 0:
            my_log.info("(excluded {} genomes from SS)".format(str(len(genomes_exclude_ss))))
        my_log.info("{} SS sequences done.".format(str(n_sequences_ss)))

    def generate_kmer_features(self):
        my_log = logging.getLogger('train:generate_kmers')
        fasta2kmers_command = utils.fasta2kmers_base(self.config, self.workingdir)
        fasta2kmers2_command = utils.fasta2kmers2_base(self.config, self.workingdir)

        use_fasta2kmers2 = int(self.config.settings["kmer_normalization"]) <= 1

        if use_fasta2kmers2:
            my_log.debug('basic command for generating kmers:\n{} -i INPUT -f OUTPUT'.format(fasta2kmers2_command))
        else:
            my_log.debug('basic command for generating kmers:\n{} INPUT >> OUTPUT'.format(fasta2kmers_command))

        tasks_by_fl = {}
        max_entries = 0
        for fl in self.config.settings["fragment_len"]:
            fl_tasks = []
            p = "{d}{sep}sampled_fasta{sep}{fl}".format(d=self.config.settings["project_dir"], fl=str(fl), sep=os.path.sep)
            outfile = "{d}{sep}train_data{sep}{fl}.sl".format(d=self.config.settings["project_dir"], fl=str(fl), sep=os.path.sep)

            files = os.listdir(p)
            for f in files:
                fastafile = os.path.join(p, f)
                if use_fasta2kmers2:
                    file_cmd = "{cmd} -i {fasta} -f {out}".format(cmd=fasta2kmers2_command, out=outfile, fasta=fastafile)
                else:
                    file_cmd = "{cmd} {fasta} >> {out}".format(cmd=fasta2kmers_command, out=outfile, fasta=fastafile)

                if self.config.settings["processors"] == 1:
                    s = os.system('{}'.format(file_cmd))
                    if s != 0:
                        my_log.critical("problem with generating kmers..:\n{}".format(file_cmd))
                        sys.exit(1)
                else:
                    fl_tasks.append(parallel.TaskCmd(file_cmd)) # attention, they all write to the same file

            if self.config.settings["processors"] > 1:
                tasks_by_fl[fl] = fl_tasks
                if len(fl_tasks) > max_entries:
                    max_entries = len(fl_tasks)

        if self.config.settings["processors"] > 1:
            # now run commands in parallel but only, if they do not write to the same file
            # this means you can only run len(fragment_len) commands in parallel
            # instead: concat files and use generate_kmer_features_concat()
            for i in range(max_entries):
                l = []
                for k in tasks_by_fl.keys():
                    try:
                        l.append(tasks_by_fl[k][i])
                    except IndexError:
                        pass
                if parallel.reportFailedCmd(parallel.runCmdParallel(l, self.config.settings["processors"])) is not None:  # Ivan change
                    sys.exit(-1)  # Ivan change

    def generate_kmer_features_concat(self):
        my_log = logging.getLogger('train:generate_kmers')
        fasta2kmers_command = utils.fasta2kmers_base(self.config, self.workingdir)
        fasta2kmers2_command = utils.fasta2kmers2_base(self.config, self.workingdir)

        use_fasta2kmers2 = int(self.config.settings["kmer_normalization"]) <= 1

        if use_fasta2kmers2:
            my_log.debug('basic command for generating kmers:\n{} -i INPUT -f OUTPUT'.format(fasta2kmers2_command))
        else:
            my_log.debug('basic command for generating kmers:\n{} INPUT >> OUTPUT'.format(fasta2kmers_command))

        tasks = []
        for fl in self.config.settings["fragment_len"]:
            p = "{d}{sep}sampled_fasta{sep}{fl}".format(d=self.config.settings["project_dir"], fl=str(fl), sep=os.path.sep)
            combined_fasta = os.path.join(p, "{}.all.fna".format(fl))
            outfile = "{d}{sep}train_data{sep}{fl}.sl".format(d=self.config.settings["project_dir"], fl=str(fl), sep=os.path.sep)
            os.system("cat {dir}{sep}*.fna > {dir}{sep}{fl}.all.tmp".format(dir=p, fl=fl, sep=os.path.sep))
            os.system("rm {dir}{sep}*.fna".format(dir=p, sep=os.path.sep))
            os.system("mv {p}{sep}{fl}.all.tmp {combined}".format(p=p, fl=fl, combined=combined_fasta, sep=os.path.sep))

            # in-efficent:
            # files = os.listdir(p)
            # s = os.system("mv {file0} {combined}".format(file0=os.path.join(p, files[0]), combined=combined_fasta))
            # if s != 0:
            #     sys.stderr.write("problem with moving file {}\n".format(os.path.join(p,files[0])))
            #     sys.exit(1)
            #
            # for f in files[1:]:
            #     s = os.system("cat {combined} {f} >> {combined}.tmp".format(dir=p, combined=combined_fasta, f=os.path.join(p, f)))
            #     if s != 0:
            #         sys.stderr.write("Problem with concatenating files in sampled_fasta/{fl}\n".format(fl=fl))
            #     os.system("mv {c}.tmp {c}".format(c=combined_fasta))
            #     os.remove(os.path.join(p, f))

            if use_fasta2kmers2:
                command = "{cmd} -i {combined} -f {out}".format(cmd=fasta2kmers2_command, combined=combined_fasta, out=outfile)
            else:
                command = "{cmd} {combined} >> {out}".format(cmd=fasta2kmers_command, combined=combined_fasta, out=outfile)

            if self.config.settings["processors"] > 1:
                tasks.append(parallel.TaskCmd(command))
            else:
                s = os.system('{}'.format(command))
                if s != 0:
                    my_log.critical("problem with generating kmers..:\n{}".format(command))
                    sys.exit(1)

        if self.config.settings["processors"] > 1:
            if parallel.reportFailedCmd(parallel.runCmdParallel(tasks, self.config.settings["processors"])) is not None:  # Ivan change
                sys.exit(-1)  # Ivan change

    def build_models(self):
        my_log = logging.getLogger('train:build_models')
        # now as the training data is ready get the models
        # if no grid was given then just build models
        # kernel options
        kernel_opt = "-t {t} -g {g} -d {d} -s {s}".format(t=str(self.config.settings["kernel"]),
                                                           g=str(self.config.settings["kernel_rbf_gamma"]),
                                                           d=str(self.config.settings["kernel_polynomial_degree"]),
                                                           s=str(self.config.settings["kernel_polynomial_s"]))
        loss_opt = "-l {l} --L {L}".format(l=str(self.loss_function),
                                            L=str(self.config.settings["loss_action"]))
        other_opt = "--z {z} --v {v} --t {t}".format(z=str(self.z_standardization),
                                                      v=str(self.misc_nodes),
                                                      t=self.tree_file)
        learn_command = "{bin} {kernel} {loss} {other} " \
                        "-v 1 -o 2".format(bin=utils.path_to_binary(self.workingdir, "svm_phylo_learn"),
                                           kernel=kernel_opt,
                                           loss=loss_opt,
                                           other=other_opt)
        cv_command = "{bin} {kernel} {loss} {other} " \
                     "-x 3 -v 1 -o 2 --r 1 --S 1".format(bin=utils.path_to_binary(self.workingdir, "svm_phylo_cv"),
                                                         kernel=kernel_opt,
                                                         loss=loss_opt,
                                                         other=other_opt)
        if self.config.settings["balance_classes"]:
            learn_command = "{} --c 1".format(learn_command)
            cv_command = "{} --c 1".format(cv_command)

        my_log.debug('basic crossvalidation command:\n{} -c CVAL KMER_FILE'.format(cv_command))
        my_log.debug('basic learning command:\n{} -c CVAL KMER_FILE MODEL_FILE'.format(learn_command))

        tasks = []  # only needed when running in parallel

        if len(self.config.settings["c_grid"]) == 1:
            c_val = self.config.settings["c_grid"][0]
            for fl in self.config.settings["fragment_len"]:
                learn_command_final = "{cmd} -c {cval} " \
                                      "{p}{sep}train_data{sep}{fl}.sl " \
                                      "{p}{sep}models{sep}{fl}_c{cval}.svm".format(cmd=learn_command,
                                                                                   sep=os.path.sep,
                                                                                   cval=c_val,
                                                                                   p=self.config.settings["project_dir"],
                                                                                   fl=fl)
                if self.config.settings["processors"] == 1:
                    my_log.info("build {fl} length model with c={cval}".format(fl=fl, cval=c_val))
                    s = os.system(learn_command_final)
                    if s != 0:
                        my_log.critical("something went wrong with building the model:\n{}".format(learn_command_final))
                        sys.exit(1)
                else:
                    tasks.append(parallel.TaskCmd(learn_command_final))

        else:
            # crossvalidation
            for fl in self.config.settings["fragment_len"]:
                my_log.info("Cross-validatong {} length model.".format(fl))
                cv_loss = []
                cv_zeroone = []
                for c_val in self.config.settings["c_grid"]:
                    my_log.debug("c={}".format(c_val))
                    cv_command_final = "{cmd} -c {cval} " \
                                       "{p}{sep}train_data{sep}{fl}.sl".format(cmd=cv_command,
                                                                               sep=os.path.sep,
                                                                               cval=c_val,
                                                                               p=self.config.settings["project_dir"],
                                                                               fl=fl)
                    fr = sys.stdin
                    os.system(cv_command_final)
                    lines = fr.readlines()
                    fr.close()
                    floating_point = re.compile(r'\d+\.\d+')
                    for line in lines:
                        if "Average loss in cross-validation" in line:
                            loss = floating_point.findall(line)
                            try:
                                loss = float(loss[0])
                                cv_loss.append(loss)
                            except IndexError:
                                continue
                        if "one-error in cross-validation" in line:
                            loss = floating_point.findall(line)
                            try:
                                loss = float(loss[0])
                                cv_zeroone.append(loss)
                            except IndexError:
                                continue

                if len(cv_loss) != len(self.config.settings["c_grid"]):
                    my_log.critical("Error, something went wrong with cross-validation "
                                    "of {} length fragment model. Quitting".format(fl))   # exit or continue?
                    sys.exit(1)
                my_log.debug("C grid: " + utils.any_list_to_string(self.config.settings["c_grid"]))
                my_log.debug("CV loss: " + utils.any_list_to_string(cv_loss))
                my_log.debug("CV 0-1: " + utils.any_list_to_string(cv_zeroone))

                loss_min = min(cv_loss)
                i = cv_loss.index(loss_min)
                # build model for minimum loss
                c_val = self.config.settings["c_grid"][i]
                learn_command_final = "{cmd} -c {cval} " \
                                      "{p}{sep}train_data{sep}{fl}.sl " \
                                      "{p}{sep}models{sep}{fl}_c{cval}.svm".format(cmd=learn_command,
                                                                                   sep=os.path.sep,
                                                                                   cval=c_val,
                                                                                   p=self.config.settings["project_dir"],
                                                                                   fl=fl)
                my_log.info("build {fl} length model with c={cval} and CV-loss={loss}".format(fl=fl,
                                                                                                cval=c_val,
                                                                                                loss=loss_min))
                if self.config.settings["processors"] == 1:
                    s = os.system("{}".format(learn_command_final))
                    if s != 0:
                        my_log.critical("something went wrong with building the model:\n{}".format(learn_command_final))
                        sys.exit(1)
                else:
                    tasks.append(parallel.TaskCmd(learn_command_final))

        if self.config.settings["processors"] > 1:
            my_log.info("building models in parallel...")
            if parallel.reportFailedCmd(parallel.runCmdParallel(tasks, self.config.settings["processors"])) is not None:  # Ivan change
                sys.exit(-1)  # Ivan change

    def get_mapped_organism(self, organism):
        """
            get the node on the tree, where the organism can be mapped to
            return None if it could not be mapped at all
        """
        mapped = None
        if str(organism) in self.nodes:
            mapped = organism
        if mapped is None:
            for rank in self.config.settings["taxonomy_ranks"]:
                parent = str(self.sqlite_taxonomy.parentAtRank(organism, rank))
                if parent in self.nodes and parent != "":
                    mapped = parent
                    break
        return mapped

    def create_folderstructure(self):
        """
            creates the folder structure for internal project directory
            project_dir/labels/
            project_dir/train_data/
            project_dir/models/
            project_dir/sampled_fasta/fragment_length_1
            ...
            project_dir/sampled_fasta/fragment_length_n
        """

        directories = ["labels", "sampled_fasta", "train_data", "models"]
        for directory in directories:
            d = os.path.join(self.config.settings["project_dir"], directory)
            os.mkdir(d)

        d = os.path.join(self.config.settings["project_dir"], "sampled_fasta")
        for fl in self.config.settings["fragment_len"]:
            subdir = os.path.join(d, str(fl))
            os.mkdir(subdir)

        #                 TAXONOMY OPTIONS                  #

    def check_db(self):
        # check if the sqlite database exists otherwise create one
        # assuming a constant name for the sqlite database
        my_log = logging.getLogger('train:database_check')
        if not os.path.isfile(self.config.ncbi_tax_db):
            my_log.info("\nThe SQLite database doesn't exist we will create one ...")
            scr = os.path.join(self.workingdir, "ncbitax2sqlite.py")
            tax_call = "python {scr} {taxdir} {taxdb}".format(scr=scr,
                                                              taxdir=self.config.settings["ncbi_tax_dir"],
                                                              taxdb=self.config.ncbi_tax_db)
            if self.yes:
                tax_call = "{} -y".format(tax_call)
            my_log.debug(tax_call)
            success = os.system(tax_call)
            if success != 0:
                my_log.critical("Error in creating the SQLite database.")
                sys.exit(1)
            my_log.debug("The NCBI database did not exist - created one: {}".format(self.config.ncbi_tax_db))
        self.sqlite_taxonomy = taxonomy_ncbi.TaxonomyNcbi(self.config.ncbi_tax_db)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="configuration file", action='store', required=True)
    parser.add_argument('-y', help="automatically set answers to 'yes'", action='store_true', default=False)
    parser.add_argument("-fb", help="backup directory from a previous PPS run, "
                                    "will try to restart from within the pipeline", action='store', default=None)
    parser.add_argument('-b', help="write a backup to the project directory", action='store_true', default=False)
    args = parser.parse_args()

    logmethods.initialise('train.log')
    cfg = config.Config(args.c, logged=True)
    obj = Train(cfg, args.y, args.b, logged=True)

    if args.fb is None:
        obj.main_processing()
    else:
        obj.from_backup(args.fb)


