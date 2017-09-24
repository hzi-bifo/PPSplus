#!/usr/bin/env python
"""
Author: Jessika Fiedler (jessika.fiedler@hhu.de)
Copyright (c) Jessika Fiedler 2015

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
    Used to store and load the pipeline from different points
        - keep track of written files
        - store the current values of each variable
    Call function write_backup() to write the complete status to a specified directory
        -> creates the following files:
            "files.dmp"             contains a Filestatus-object for all successfully written files
            "failed_files.dmp"      contains Filestatus-objects for files that have not been written successfully
                                        (does not exist if all files have been written correctly)
            "{variable}.dmp"         for each stored variable with the current value
"""


import pickle
import os
import time
import sys
import logging
import logmethods


class Status:

    def __init__(self, logged=False):
        self.verbose = False

        self.variables = {}
        self.files = {}     # path -> Filestatus

        # only needed when a status is loaded from backup
        self.failed_files = {}

        if not logged:
            # this will write a file to cwd
            logmethods.initialise('status.log')

        self.log = logging.getLogger('status')

    def add_written_file(self, f):
        """
           whenever a file is written, store its path
            used to check whether all files have been written successfully
            and to check if all written files are still available, when loading an object from this status
        """
        self.log.debug('added file:\t{}'.format(f))
        self.files[f] = Filestatus(f)

    def succesfully_written(self, f):
        self.log.debug('closed file:\t{}'.format(f))
        self.files[f].closed()

    def change_variable(self, obj, variable_name):
        """
            should be called for all necessary and modified variables whenever the program reaches a critical part
        """
        self.log.debug('changed variable:\t{}'.format(variable_name))
        self.variables[variable_name] = obj

    def write_backup(self, backupdir, rm_old=True):
        self.log.debug('writing backup to:\n{}'.format(backupdir))

        if rm_old:
            try:
                self.log.info('overwriting last backup')
                shutil.rmtree(backupdir)
                os.mkdir(backupdir)
            except:
                self.log.warning('failed to delete backup from last run')

        not_closed = []
        with open(os.path.join(backupdir, "files.dmp"), 'wb') as output:
            for f in self.files:
                if not self.files[f].successfully_written():
                    self.log.debug('file has not been successfully written:\n{}'.format(f))
                    not_closed.append(f)
                else:
                    pickle.dump(self.files[f], output, pickle.HIGHEST_PROTOCOL)

        if len(not_closed) > 0:
            self.log.debug('there are files that have not been successfully written. Backup to failed_files.dmp')
            with open(os.path.join(backupdir, "failed_files.dmp"), 'wb') as output:
                for f in not_closed:
                    pickle.dump(self.files[f], output, pickle.HIGHEST_PROTOCOL)

        for v in self.variables:
            with open(os.path.join(backupdir, "{}.dmp".format(v)), 'wb') as output:
                pickle.dump(self.variables[v], output, pickle.HIGHEST_PROTOCOL)

    def load_backup(self, backupdir):
        """
            loads a backup and checks its consistency
                - the written files should still have the same modification and creation times
                - files should not be empty
                - if a file has not been written successfully, it is stored to failed_files
        """
        files = os.listdir(backupdir)
        self.log.debug('loading backup from:\n{}'.format(backupdir))
        if "failed_files.dmp" in files:
            files.remove("failed_files.dmp")
            with open(os.path.join(backupdir, 'failed_files.dmp'), 'rb') as input:
                while True:
                    try:
                        f = pickle.load(input)
                    except:
                        break
                    self.failed_files.append(f)

        if "files.dmp" in files:
            files.remove('files.dmp')
            with open(os.path.join(backupdir, 'files.dmp'), 'rb') as input:
                while True:
                    try:
                        f = pickle.load(input)
                    except:
                        break

                    m_time = time.ctime(os.path.getmtime(f.path))
                    c_time = time.ctime(os.path.getctime(f.path))
                    if not f.same_time(m_time, c_time):
                        self.log.critical("The file \n\t{}\n has been modified since the last run\n"
                                          "Could not load the pipeline correctly.".format(f.path))
                        sys.exit(1)

                    f_obj = open(f.path)
                    if len(f_obj.readlines()) == 0:
                        self.failed_files.append(f)
                        f_obj.close()
                        continue
                    f_obj.close()

                    # file is not empty and has same modification and creation time
                    self.files[f.path] = f

        # now load the variables
        for basename in files:
            variable_name = basename.split(".")[0]
            with open(os.path.join(backupdir, basename), 'rb') as input:
                self.variables[variable_name] = pickle.load(input)


class Filestatus:

    def __init__(self, f):
        self.open = True
        self.last_modified = time.ctime(os.path.getmtime(f))
        self.creation_time = time.ctime(os.path.getctime(f))
        self.path = f

    def closed(self):
        self.open = False

    def successfully_written(self):
        # if a file has been successfully written, it is closed ( = not open)
        return not self.open

    def same_time(self, mtime, ctime):
        if mtime == self.last_modified and ctime == self.creation_time:
            return True
        else:
            return False

import unittest
import tempfile
import shutil


class StatusTesting(unittest.TestCase):

    def test_written_file(self):
        directory = tempfile.mkdtemp()
        backupdir = os.path.join(directory, "backup")

        s = Status()
        filepath = os.path.join(directory, "tempfile.txt")

        lines = []
        for i in range(20):
            lines.append(str(i))

        t_file = open(filepath, "w")
        s.add_written_file(filepath)

        t_file.writelines(lines)
        t_file.close()
        s.succesfully_written(filepath)

        os.mkdir(backupdir)
        s.write_backup(backupdir)

        s2 = Status()
        s2.load_backup(backupdir)

        self.assertListEqual(s.files.keys(), s2.files.keys())
        self.assertListEqual(s.variables.keys(), s2.variables.keys())
        self.assertListEqual(s.failed_files.keys(), s2.failed_files.keys())

        os.utime(filepath, (1330712280, 1330712292))
        s1 = Status()
        with self.assertRaises(SystemExit) as cm:
            s1.load_backup(backupdir)
        self.assertEqual(cm.exception.code, 1)

        shutil.rmtree(backupdir)

    def test_written_files(self):
        directory = tempfile.mkdtemp()
        backupdir = os.path.join(directory, "backup")

        s = Status()
        s2 = Status()

        for j in range(10):
            filepath = os.path.join(directory, "tempfile{}.txt".format(str(j)))
            lines = []
            for i in range(20):
                lines.append(str(i))

            t_file = open(filepath, "w")
            s.add_written_file(filepath)

            t_file.writelines(lines)
            t_file.close()
            s.succesfully_written(filepath)

        os.mkdir(backupdir)
        s.write_backup(backupdir)

        s2.load_backup(backupdir)

        self.assertListEqual(s.files.keys(), s2.files.keys())
        self.assertListEqual(s.variables.keys(), s2.variables.keys())
        self.assertListEqual(s.failed_files.keys(), s2.failed_files.keys())

        os.utime(os.path.join(directory, "tempfile3.txt".format(str(j))), (1330712280, 1330712292))
        s1 = Status()
        with self.assertRaises(SystemExit) as cm:
            s1.load_backup(backupdir)
        self.assertEqual(cm.exception.code, 1)

        shutil.rmtree(backupdir)

    def test_write_variables(self):
        directory = tempfile.mkdtemp()
        backupdir = os.path.join(directory, "backup")

        v = dict()
        v["v1"] = True
        v["v2"] = 10
        v["v3"] = "hallo"
        v["v4"] = [1, 2, 3, 4]
        v["v5"] = ["a", "b", "c"]

        s = Status()
        for v_name in v:
            s.change_variable(v[v_name], v_name)

        self.assertDictEqual(s.variables, v)

        os.mkdir(backupdir)
        s.write_backup(backupdir)

        s_loaded = Status()
        s_loaded.load_backup(backupdir)

        self.assertDictEqual(s.variables, s_loaded.variables)

        shutil.rmtree(backupdir)

    if __name__ == '__main__':
        unittest.main()