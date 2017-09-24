#!/usr/bin/env python

"""
    Copyright (C) 2014  Jessika Fiedler
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
    Script to build a sqlite3 database from NCBI data.
"""

import sqlite3
import argparse
import os
import time
import sys


def get_answer_timeout():
    start_time = time.time()  # this is time in seconds
    while True:
        answer = sys.stdin.readline().strip().lower()
        if answer == "y":
            return "y"
        curr_time = time.time()
        time_ellapsed = curr_time - start_time
        if time_ellapsed is 120 or answer == "n":
            return "n"


def build_database(args):
    """
            TABLE taxon
                taxon_id            autoincrementing integer
                ncbi_taxon_id       the ncbi id
                parent_taxon_id     ncbi id of the taxons parent
                node_rank           the taxons rank
            TABLE taxon_name
                taxon_id    same as in TABLE taxon
                name        the taxons name at given namespace
                name_class  e.g. 'scientific name'
    """
    if os.path.isfile(args.db):
        database_exists()
    # create a new database
    db = sqlite3.connect(args.db)
    cursor = db.cursor()
    # specify tables in the database
    taxon_table = "CREATE TABLE taxon(" \
                  "taxon_id INTEGER PRIMARY KEY AUTOINCREMENT," \
                  "ncbi_taxon_id INTEGER," \
                  "parent_taxon_id  INTEGER," \
                  "node_rank TEXT," \
                  "UNIQUE (ncbi_taxon_id));"
    cursor.execute(taxon_table)
    cursor.execute("CREATE INDEX taxparent ON taxon(parent_taxon_id);")
    name_table = "CREATE TABLE taxon_name(" \
                 "taxon_id INTEGER," \
                 "name TEXT NOT NULL," \
                 "name_class TEXT NOT NULL," \
                 "UNIQUE (taxon_id, name, name_class));"
    cursor.execute(name_table)

    sys.stdout.write("Processing nodes.dmp... \n")
    # read the ncbi dumps and populate the database
    # first load the nodes.dmp as we will need it for names.dmp
    # using parameterized queries for performance
    # examples at http://rosettacode.org/wiki/Parametrized_SQL_statement
    fr = open(os.path.join(args.dmp, "nodes.dmp"))
    for line in fr:
        line = line.strip()
        if line == "":
            continue
        values = line.split("|", 3)

        nid = values[0].strip("\t")
        pid = values[1].strip("\t")
        rank = values[2].strip("\t")
        insert_taxon = "INSERT INTO taxon (ncbi_taxon_id, parent_taxon_id, node_rank) " \
                       "VALUES({ncbi_id},{parent_id},'{rank}')".format(ncbi_id=nid, parent_id=pid, rank=rank)
        cursor.execute(insert_taxon)
    fr.close()
    sys.stdout.write("Done.\n")

    sys.stdout.write("Processing names.dmp... \n")
    fr = open(os.path.join(args.dmp, "names.dmp"))
    for line in fr:
        line = line.strip()
        if line == "":
            continue
        values = line.split("|", 4)
        ncbi_taxon_id = values[0].strip("\t")
        cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?', (ncbi_taxon_id,))
        result = cursor.fetchall()
        taxon_id = str(result[0][0])
        name = values[1].strip("\t")
        if "\'" in name:
            name = name.replace('\'', "\'\'")   # assuming '' in sql has same effect as \' in python
        name_class = values[3].strip("\t")
        insert_taxon_name = "INSERT OR IGNORE INTO taxon_name VALUES({taxid},'{name}','name_class')".format(
            taxid=taxon_id, name=name, name_class=name_class)
        cursor.execute(insert_taxon_name)

    fr.close()
    db.commit()
    db.close()
    sys.stdout.write("Done.")


def database_exists(checkold= True):
    # check if the database is valid
    if checkold:
        check_passed = False
        db = sqlite3.connect(args.db)
        cursor = db.cursor()
        try:
            cursor.execute('SELECT name FROM taxon_name T WHERE T.taxon_id=?', (1,))
            res1 = cursor.fetchall()
            cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.taxon_id=?', (1,))
            res2 = cursor.fetchall()
            if len(res1) == 1 and len(res2) == 1:
                check_passed = True
        except:
            check_passed = False
        if check_passed:
            sys.stdout.write("It passed the check and is a valid database (possibly outdated).")
        else:
            sys.stdout.write("It did not pass the check and is corrupted.")
    else:
        check_passed = True
        # it is a 'simple' database
        # Todo: implement check

    if args.y:
        ans = "y"
    else:
        sys.stdout.write("Remove it? [Y/N] (default=Y, timeout 2 minutes)")
        ans = get_answer_timeout()

    if ans is "y":
        sys.stdout.write("removing old database...")
        os.remove(args.db)
    else:
        if check_passed:
            sys.stdout.write("The database will not be touched.")
        else:
            sys.stdout.write("Please run again with valid options, abort.")
        sys.exit(1)


def build_database_simple(args):
    """
            TABLE taxon_simple
                ncbi_taxon_id       the ncbi id
                parent_taxon_id     ncbi id of the taxons parent
                rank                the taxons rank
                scientific_name     the taxons name at given namespace
        Todo: read in dict first, then create db
    """

    if os.path.isfile(args.db):
        database_exists(checkold=False)

    # create a new database
    db = sqlite3.connect(args.db)
    cursor = db.cursor()
    # specify table in the database
    taxon_table = "CREATE TABLE taxon_simple(" \
                  "ncbi_taxon_id INTEGER PRIMARY KEY," \
                  "parent_ncbi_taxon_id  INTEGER NOT NULL," \
                  "rank TEXT NOT NULL," \
                  "scientific_name TEXT);"
    cursor.execute(taxon_table)
    cursor.execute("CREATE INDEX taxon_simple_parent_index ON taxon_simple(parent_ncbi_taxon_id)")

    sys.stdout.write("Processing nodes.dmp... ")
    # read the ncbi dumps and populate the database
    # first load the nodes.dmp as we will need it for names.dmp
    fr = open(os.path.join(args.dmp, "nodes.dmp"))
    for line in fr:
        if line == "":
            continue
        values = line.split("|")[0:3]
        insert_taxon = "INSERT INTO taxon_simple (ncbi_taxon_id, parent_ncbi_taxon_id, rank) " \
                       "VALUES({nid},{pid},'{rank}')".format(nid=values[0].strip("\t"),
                                                             pid=values[1].strip("\t"),
                                                             rank=values[2].strip("\t"))
        cursor.execute(insert_taxon)
    fr.close()
    sys.stdout.write("Done.")

    sys.stdout.write("Processing names.dmp... ")
    fr = open(os.path.join(args.dmp, "names.dmp"))
    for line in fr:
        if line == "":
            continue
        values = line.strip().split("|")
        name_class = values[3].strip("\t")
        # store scientific names only
        if name_class == 'scientific name':
            ncbi_taxon_id = values[0].strip("\t")
            name = values[1].strip("\t")
            if "\'" in name:
                name = name.replace('\'', "\'\'")   # assuming '' in sql has same effect as \' in python

            # use update to insert value into existing database:
            update = "UPDATE taxon_simple " \
                     "SET scientific_name = '{name}' " \
                     "WHERE ncbi_taxon_id = {nid}".format(name=name, nid=ncbi_taxon_id)
            cursor.execute(update)
    fr.close()
    db.commit()
    db.close()
    sys.stdout.write("Done.")


def download_dumps(args):
    sys.stdout.write("\tdownloading dump files...\n")  # TODO: use python, not wget, to download and only extract required files

    data_archive = os.path.join(args.dmp, "taxdump.tar.gz")

    download_cmd = "wget -O {} ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz".format(data_archive)
    sucess = os.system(download_cmd)
    if sucess != 0:
        sys.stderr.write("Something went wrong with downloading ncbi data\n")
        sys.exit(1)
    sys.stdout.write("\tunpacking data...")
    unpack_cmd = "tar xfz {source} -C {dmp_dir}".format(source=data_archive, dmp_dir=args.dmp)
    sucess = os.system(unpack_cmd)
    if sucess != 0:
        sys.stderr.write("Something went wrong with unpacking data\n")
        sys.exit(1)
    os.remove(data_archive)
    if sucess != 0:
        sys.stdout.write("failed to remove {}\n".format(data_archive))

if __name__ == "__main__":
    '''
        Uses NCBI nodes.dmp and names.dmp to fill the database
        If not provided, the NCBI taxonomy dumps are downloaded from the NCBI FTP site.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument("-dmp", help="directory including ncbi dumpfiles", action='store', required=True)
    parser.add_argument("-db", help="filename for the SQLite database", action='store', required=True)
    parser.add_argument('-y', help="automatically set answers to 'yes'", action='store_true', default=False)
    parser.add_argument('-s', help="build a more simple variant of the database", action='store_true', default=False)
    args = parser.parse_args()

    if not os.path.isfile(os.path.join(args.dmp, "nodes.dmp")) or not os.path.isfile(os.path.join(args.dmp, "names.dmp")):
        sys.stdout.write("NCBI taxonomy dump files are not present.\n")
        if args.y:
            ans = "y"
        else:
            sys.stdout.write("Download taxonomy data from NCBI? [Y/N] (default=Y, timeout 2 minutes)\n")
            ans = get_answer_timeout()

        if ans is "y":
            download_dumps(args)
        else:
            sys.stderr.write("\nPlease provide correct directory name with NCBI taxonomy dump files.\n")
            sys.exit(1)

    if args.s:
        build_database_simple(args)
    else:
        build_database(args)

