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

About: Parser for PPS configuration file
"""

import sys
import os
import multiprocessing as mp
from taxonomy_ncbi import TAXONOMIC_RANKS
import logging
import logmethods

class Config:
    def __init__(self, config_file, values=None, logged=False):
        """
            reads all options in the configuration file and stores them into a dictionary
            checks values and presence of necessary information
            instead of passing a configuration file, you can pass the values directly as a dictionary
        """
        if not logged:
            # this will write a file to cwd
            logmethods.initialise('config.log')
        self.log = logging.getLogger('config')

        self.settings = {}

        self.config_script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
        self.tax_string = ""
        self.ncbi_tax_db = None
        self.genomes_exclude = []

        if config_file is None:
            if values is not None:
                self.settings = values
                self.tax_string = ",".join(self.settings["taxonomy_ranks"])
            else:
                self.log.critical("Please provide either values or a configuration file.")
                sys.exit(1)
        else:
            self._read_file(config_file)

        self._check_values()

    def _read_file(self, config_file):
        self.log.debug("Processing according to {}".format(config_file))
        self.log.info("reading configuration file...")

        for line in open(config_file):
            line = line.strip()
            if line is None or line.startswith("#") or line == "":
                continue
            if len(line.split(":")) == 2:
                key, value = line.split(":")
                if key != "":
                    self.settings[key.lower()] = value
            elif len(line.split(":")) > 2:  # happens when using windows paths for example
                key = line.split(":")[0]
                value = line.strip()[line.find(":")+1:]
                if key != "":
                    self.settings[key.lower()] = value
        try:
            self.tax_string = self.settings["taxonomy_ranks"]
        except KeyError:
            TAXONOMIC_RANKS.reverse()
            self.tax_string = ",".join(TAXONOMIC_RANKS)
            self.settings["taxonomy_ranks"] = self.tax_string

        self._config_store_list("fragment_len", ",", convert_to=int)
        self._config_store_list("c_grid", ",", convert_to=int)
        self._config_store_list("kmer", "-", convert_to=int)
        self._config_store_list("taxonomy_ranks", ",")
        self._config_store_list("extensions", ",", default=[])
        self._config_store_list("sample_specific_step", ",", convert_to=int, default=[])
        self._config_store_number("number_examples", float)
        self._config_store_number("n_min_genomes_generic", int)
        self._config_store_number("kernel", int)
        self._config_store_number("kernel_polynomial_degree", int)
        self._config_store_number("kernel_polynomial_s", float)
        self._config_store_number("kernel_rbf_gamma", float)
        self._config_store_number("kmer_normalization", int)
        self._config_store_number("rev_complement", int)
        self._config_store_number("rm_rev_complement", int, default=0)
        self._config_store_number("loss_action", int)
        self._config_store_number("processors", int, default=mp.cpu_count())
        self._config_store_bool("only_models")
        self._config_store_bool("balance_classes")
        self._config_store_bool("postprocessing")
        self._config_store_bool("clean_up_train")
        self._config_store_number("n_classifiers", int)

    def _check_values(self):
        # check if necessary information is provided:
        necessary = ["ncbi_tax_dir", "ncbi_processed_dir", "project_dir", "fragment_len",
                     "rev_complement", "kmer_normalization", "loss_action"]
        for k in necessary:
            if not self._key_is_present(k):
                self.log.critical("please define option {} in the configuration file.".format(k))
                sys.exit(1)
        # check if directories are valid
        for d in ["ncbi_processed_dir", "project_dir", "ncbi_tax_dir"]:
            if not os.path.isdir(self.settings[d]):
                self.log.critical("Invalid {} entry in the configuration file".format(d))
                sys.exit(1)
        # valid loss action?
        if self.settings["loss_action"] != 0 and self.settings["loss_action"] != 1:
            self.log.critical("Invalid loss action. Expected 0 or 1. Quitting..")
            sys.exit(1)
        # valid kmer_normalization
        if self.settings["kmer_normalization"] < 0 or self.settings["kmer_normalization"] > 3:
            self.log.critical("Invalid kmer_normalization. Expected 0,1,2 or 3. Quitting..")
            sys.exit(1)
        # if a tree file is provided, check if it is a valid file
        if self._key_is_present("tree_file"):
            if not os.path.isfile(self.settings["tree_file"]):
                self.log.critical("Invalid TREE_FILE entry in the configuration file")
                sys.exit(1)
        else:
            # if no tree file is provided check if a valid number of generic genomes to build a tree from is given
            if not self._key_is_present("n_min_genomes_generic") or self.settings["n_min_genomes_generic"] <= 0:
                self.log.critical("Invalid N_MIN_GENOMES_GENERIC entry, must be a positive integer")
                sys.exit(1)
            else:
                self.log.info("Tree file not provided, will build a generic model with {} min clade genomes"
                              .format(str(self.settings["n_min_genomes_generic"])))
                self.settings["tree_file"] = None
        # check if all necessary lists have entries
        for l in ["fragment_len", "taxonomy_ranks", "c_grid"]:
            if len(self.settings[l]) == 0:
                self.log.critical("Invalid {} entry in the configuration file".format(l))
                sys.exit(1)
        # if a directory with sample specific data is given, check if it is a valid directory
        if self._key_is_present("sample_specific_dir") and not os.path.isdir(self.settings["sample_specific_dir"]):
            self.log.critical("Invalid SAMPLE_SPECIFIC_DIR entry in the configuration file")
            sys.exit(1)

            # for each fragment length, there must be a sample specific step
            # if there is only one value given as sample specific step,
            # you store it in a list with the same length as fragment length
            #    fragment lengths           [a,b,c]
            #    sample specific steps      [d,e,f] or [d,d,d]
            nr_of_sample_specific_steps = len(self.settings["sample_specific_step"])
            if nr_of_sample_specific_steps == 1:
                self.settings["sample_specific_step"] *= len(self.settings["fragment_len"])
                # e.g. 3*[c] = [c,c,c]
            elif nr_of_sample_specific_steps > 1 and (nr_of_sample_specific_steps != len(self.settings["fragment_len"])):
                self.log.critical("Not same number of SAMPLE_SPECIFIC_STEP for each FRAGMENT_LEN")
                sys.exit(1)

        # if there is a list of genomes which should be excluded from training,
        # check if the file is valid and store them in a list
        if self._key_is_present("genomes_exclude"):
            if not os.path.isfile(self.settings["genomes_exclude"]):
                self.log.critical("Invalid GENOMES_EXCLUDE entry in the configuration file")
                sys.exit(1)
            else:
                for line in open(self.settings["genomes_exclude"]):
                    line = line.strip()
                    if line != "":
                        self.genomes_exclude.append(line)
        # store the complete path to the ncbi database either use the user defined name or set it to default
        if self._key_is_present("ncbi_sql"):
            self.ncbi_tax_db = os.path.join(self.settings["ncbi_tax_dir"], self.settings["ncbi_sql"])
        else:
            self.ncbi_tax_db = os.path.join(self.settings["ncbi_tax_dir"], "ncbitax_sqlite.db")

    def _config_store_bool(self, key, default=False):
        if key not in self.settings:
            # print "can not find option {} in the configuration file, will set it to False".format(key)
            self.settings[key] = default
        elif self.settings[key].lower() == "true":
            self.settings[key] = True
        elif self.settings[key].lower() == "false":
            self.settings[key] = False
        else:
            self.log.critical("expected true or false for key: {k}".format(key))
            sys.exit(1)

    def _config_store_list(self, key, sep, convert_to=None, default=None):
        if key not in self.settings:
            if default is not None:
                self.settings[key] = default
            else:
                self.log.critical("missing option {} in the configuration file".format(key))
                sys.exit(1)
        else:
            self.settings[key] = self.settings[key].split(sep)
        if convert_to is not None:
            self.settings[key] = map(lambda x: convert_to(x), self.settings[key])

    def _config_store_number(self, key, t, default=None):
        if not self._key_is_present(key):
            if default is not None:
                self.settings[key] = default
            else:
                self.log.critical("missing option {} in the configuration file".format(key))
                sys.exit(1)
        elif key not in self.settings and default is not None:
            self.settings[key] = default
        else:
            self.settings[key] = t(self.settings[key])

    def _key_is_present(self, key):
        if key in self.settings and not self.settings[key] == "":
            return True
        return False

    def get_value(self, key):
        if key not in self.settings:
            return None
        else:
            return self.setttings[key]

    def get_tax_string(self):
        return self.tax_string


def write_config_file(filepath, settings):
    list_with_comma = ["fragment_len", "taxonomy_ranks", "c_grid", "sample_specific_step", "extensions"]

    fw = open(filepath, "w")
    for k in settings.keys():
        value = settings[k]
        if k == "kmer":
            fw.write("{k}:{v1}-{v2}\n".format(k=k, v1=value[0], v2=value[1]))
        elif k in list_with_comma:
            s = ""
            for e in settings[k]:
                s += str(e)
                s += ","
            s = s[:-1]
            fw.write("{k}:{v}\n".format(k=k, v=s))
        else:
            fw.write("{k}:{v}\n".format(k=k, v=str(settings[k])))
    fw.close()


#       Testing     
import unittest
import tempfile
import shutil


class ConfigTesting(unittest.TestCase):

    def setUp(self):
        self.directory = tempfile.mkdtemp()
        self.configfile = os.path.join(self.directory, "test.config")
        fw = open(os.path.join(self.directory, "tree.txt"), "w")
        fw.close()

        for d in ["processed_dir", "taxonomy", "project_dir"]:
            os.mkdir(os.path.join(self.directory, d))

        # all the settings that MUST be present
        self.settings = {"ncbi_processed_dir": os.path.join(self.directory, "processed_dir"),
                    "ncbi_tax_dir": os.path.join(self.directory, "taxonomy"),
                    "project_dir": os.path.join(self.directory, "project_dir"),
                    "tree_file": os.path.join(self.directory, "tree.txt"),
                    "fragment_len": [1000, 3000, 5000, 10000, 15000, 50000],
                    "c_grid": [1, 10],
                    "rev_complement": 1,
                    "kmer_normalization": 1,
                    "number_examples": 10000,
                    "loss_action": 1,
                    "kmer": [1, 4],
                    "kernel": 0,
                    "kernel_polynomial_degree": 2,
                    "kernel_polynomial_s": 1.0,
                    "n_min_genomes_generic": 3,
                    "kernel_rbf_gamma": 1.0,
                    "kmer_normalization": 1,
                    "n_classifiers": 3}

    def test_read_file(self):
        # write a config file with all necessary options
        write_config_file(self.configfile, self.settings)

        strings = ["ncbi_processed_dir", "ncbi_tax_dir", "project_dir", "tree_file"]
        listings = ["fragment_len", "c_grid"]
        integers = ["rev_complement", "kmer_normalization", "number_examples", "loss_action",
                    "kernel", "kernel_polynomial_degree", "kernel_polynomial_s",
                    "n_min_genomes_generic", "kernel_rbf_gamma", "kmer_normalization", "n_classifiers"]
        # read from the file, the settings should be the same
        c = Config(self.configfile)
        read_settings = c.settings
        for k in strings:
            self.assertEqual(self.settings[k], read_settings[k])
        for k in integers:
            self.assertEqual(self.settings[k], read_settings[k])
        for k in listings:
            self.assertListEqual(self.settings[k], read_settings[k])

        shutil.rmtree(self.directory)

    def test_missing_options(self):
        # leaving out one of the basic settings should always lead to a system exit
        for k in self.settings.keys():
            settings_copy = self.settings.copy()
            settings_copy.remove(k)
            write_config_file(self.configfile, self.settings)
            with self.assertRaises(SystemExit) as cm:
                Config(self.configfile)
            self.assertEqual(cm.exception.code, 1)
            os.remove(self.configfile)
        shutil.rmtree(self.directory)

    if __name__ == '__main__':
        unittest.main()
