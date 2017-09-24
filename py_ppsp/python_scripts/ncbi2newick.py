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
outputs a tree in newick format using ncbi taxonomy
input is the node names in a file with one node per line
"""

import sys
import argparse
import taxonomy_ncbi
import os
import re
import logging
import logmethods

class Node():
    def __init__(self, label, nodelist=None):
        """
            @type nodelist: list of Node
            @type label: str
        """
        assert nodelist is None or len(nodelist) > 0
        self.label = label
        self.nodeList = nodelist
        self.parent = None

    def is_leaf(self):
        if self.nodeList is None:
            return True
        return False

    def get_child_list(self):
        assert self.nodeList is not None
        return self.nodeList

    def add_child(self, node):
        if self.nodeList is None:
            self.nodeList = [node]
        else:
            self.nodeList.append(node)


class Ncbi2Newick:
    def __init__(self, db, logged=False):
        assert os.path.isfile(db)
        self.tax = taxonomy_ncbi.TaxonomyNcbi(db)
        self.tree = None
        self.leaves = None
        self.exc = None
        if not logged:
            # this will write a file to cwd
            logmethods.initialise('Ncbi2Newick.log')
        self.log = logging.getLogger('ncbi2newick')


    def close(self):
        self.tax.close()

    def leaves_from_list(self, cladelist):
        assert os.path.isfile(cladelist)
        self.leaves = nodes_from_list(cladelist, self.log)

    def except_from_list(self, exceptlist):
        assert os.path.isfile(exceptlist)
        self.exc = nodes_from_list(exceptlist, self.log)

    def tree_from_nodes(self, nodelist, allranks=False, exc_list=None):
        assert os.path.isfile(nodelist)
        assert exc_list is None or os.path.isfile(exc_list)
        assert allranks is False or allranks is True

        self.leaves_from_list(nodelist)

        if exc_list is not None:
            self.except_from_list(exc_list)
        # get a tree for given leaves
        # allowed_ranks are allowed ranks in the taxonomy (should be bottom up,
        # e.g. ["genus","family","order","class","phylum","superkingdom"])
        # allranks=T gives equal branch lengths
        # add_nodes adds given nodes to the path irrespective of their rank
        if self.leaves is None or len(self.leaves) == 0:
            self.log.critical('no nodes to proecess\n')
            sys.exit(1)
        # create a set including all tax_ids in the topology on each rank
        dummy_node = 0
        dummy_parents = {}  # if a node does not have a specified taxonomic parent at a rank, add a dummy node

        all_nodes = set()   # all unique nodes at each rank of the topology
        for tax_id in self.leaves:
            all_nodes.add(tax_id)
            last = tax_id       # needed to store dummy parents
            for r in taxonomy_ncbi.TAXONOMIC_RANKS:
                p = self.tax.parentAtRank(tax_id, r)
                if p is None:
                    if allranks:
                        dummy_node -= 1
                        all_nodes.add(str(dummy_node))
                        dummy_parents[last] = str(dummy_node)
                        last = dummy_node
                else:
                    all_nodes.add(p)
                    last = p

        nodes = {}
        for tax_id in all_nodes:
            if tax_id == "" or tax_id == 1:
                continue
            this_node = Node(label=str(tax_id))
            parent_id = self.tax.getParentNcbid(tax_id)
            if parent_id is None:
                parent_id = dummy_parents[tax_id]
            nodes[tax_id] = this_node
            this_node.parent = parent_id

        # add the root node:
        nodes[1] = Node(label="1")
        nodes[1].parent = 1
        root_node = None
        for node_id, this_node in nodes.iteritems():
            if node_id == this_node.parent:
                root_node = this_node
            else:
                parent_node = nodes[this_node.parent]
                parent_node.add_child(this_node)

        if root_node is None:
            self.log.critical('Did not find the root of the tree.. quitting')
            sys.exit(1)

        self.tree = root_node

    def tree_to_file(self, tree_file):
        newickstring = "{};".format(get_newick(self.tree))
        self.log.debug(newickstring)
        fw = open(tree_file, "w")
        fw.write(newickstring)
        fw.close()


def get_newick(node):
    """
        Get a tree in a newick format, call with a root of the tree.
        @type node: Node
    """
    if node.is_leaf():
        return node.label
    else:
        child_newick_list = []
        for child in node.get_child_list():
            child_newick_list.append(get_newick(child))
        return '({child_newick}){node_label}'.format(child_newick=','.join(child_newick_list), node_label=node.label)


def get_nodes_from_newick(treefile):
    """
        returns all node names represented in a newick tree string
    """
    assert os.path.isfile(treefile)

    fr = open(treefile)
    newick_string = fr.readline().rstrip()
    fr.close()
    pattern = re.compile(r"[0-9]+")
    node_names = set(pattern.findall(newick_string))
    return node_names


def nodes_from_list(nodelist, log):
    assert os.path.isfile(nodelist)

    fr = open(nodelist)
    nodes = set()   # unique nodes only
    for line in fr:
        if line != "":
            nodes.add(line.rstrip())
    fr.close()
    if len(nodes) < 1:
        log.critical("No nodes to process.. Quitting..")
        sys.exit(1)
    return nodes

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-db", help="full path to NCBI database", action='store', required=True)
    parser.add_argument("-cladelist", help="file including clades (one per line)", action='store', required=True)
    parser.add_argument("-out", help="path for the output", action='store', required=True)
    parser.add_argument("-exc", help="file including clades that must be there (one per line)",
                        action='store', default=None)

    args = parser.parse_args()

    obj = Ncbi2Newick(args.db)
    obj.tree_from_nodes(args.cladelist, allranks=False, exc_list=args.exc)
    obj.tree_to_file(args.out)
    obj.close()
