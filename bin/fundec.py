#!/usr/bin/env python
###############################################################################
#
# fundec.py <tree> <database_file> <%id cutoff>
# I should probably include some regular expression recognition to make the
# counts of each annotation more robust. I imagine that as is this
# code will only cluster the major branches. I should really set the
# %ID cutoff high to avoid over clustering. I don't want that to come
# back to bite.
#
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
 
__author__ = "Joel Boyd"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
 
###############################################################################

# System imports
import sys
import argparse
import logging
import os
import shutil
import subprocess

from skbio.tree import TreeNode
from collections import Counter

# Local imports
from reroot import Rerooter
from tree_clust import Cluster
from read_database import ReadDatabase


debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################
############################### - Exceptions - ################################

class BadTreeFileException(Exception):
    pass

###############################################################################
################################ - Classes - ##################################

class FunDec:
    
    def __init__(self, db_name):
        self.db_name = db_name
        
    # ------------------------------- Parsing ------------------------------- #
    def _find_format(self, database_file):
        '''
        Quickly determine if the database file is in tsv or csv format
        
        Parameters
        ----------
        database_file: str
            path to database file
            
        Returns
        -------
        The character deliminator, either ',' or '\t'
        '''
        first_line=open(database_file).readline()
        if '\t' in first_line:
            delim = '\t'
        elif ',' in first_line: # Not even remotely robust.
            delim = ','
        return delim
        
        
        
    def _open_tree(self, tree_path):
        '''
        Open a tree file, determine what decorations are already present. Strip
        Unwanted decoration
        
        Parameters
        ----------
        tree_path: str
            Path to a file containing a phylogenetic tree, in Newick format.
            
        Returns
        -------
        skbio TreeNode object
        '''
        tree_obj=TreeNode.read(open(tree_path))
        bootstrapped = True
        for node in tree_obj.non_tips():
            if node.name:
                try:
                    float(node.name)
                except:
                    logging.debug("Tree is decorated already. Stripping all \
    previous decoration from the tree.")
                    bootstrapped = False
                    tree_obj = self._strip_tree(tree_obj)
                    break
            else:
                if bootstrapped:
                    logging.warning("This tree doesn't appear correctly \
formatted or there is information missing. No boostrap value or decoration \
found for bare node. ")
                    bootstrapped = False
        if bootstrapped:
            logging.debug("Tree is bootstrap or has confidence values \
assigned to the nodes.")
        return tree_obj

    def _strip_tree(self, tree):
        '''
        Remove all current node annotations on the tree (including bootstraps
        atm, sorry about that.)
        
        Parameters
        ----------
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode        
        Returns
        -------
        'naked' tree
        '''
        for node in tree.traverse():
            if node.is_tip:
                continue
            else:
                node.name = ''
        return tree

    def main(self, tree_path, output_path, database_path, 
             percentile, context_distance):
        
        # Setup
        tree = self._open_tree(tree_path)
        c = Cluster()
        Rerooter().reroot(tree)    
            
        if percentile:
            for index, cutoff in enumerate(percentile):
                logging.info("Identifying depth clusters in tree using %ith percentile\
 cutoff" % cutoff)
                tree=c.depth_first_cluster(tree, cutoff, index)
        
        if database_path:
            delim = self._find_format(database_path)
            database = ReadDatabase(open(database_path), delim)
            
            
            for column in database.columns():
                field = database.field(column)
                logging.info("Decorating using %s" % (field.name))
                if field.type == ReadDatabase._TYPE_NUMERIC:
                    tree=c.numeric_decoration(tree, 
                                              field.entries,
                                              field.name)
                elif field.type == ReadDatabase._TYPE_STRING:
                    tree=c.string_decoration(tree, 
                                             field.entries, 
                                             field.name)
                elif field.type == ReadDatabase._TYPE_TAXONOMY:
                    tree=c.taxonomy_decoration(tree, field.entries)
                
        logging.info("Finished clustering! Writing tree to file: %s" \
                                                            % (output_path))
        tree.write(output_path, format = "newick")       


############################### - Functions - #################################

def check_args(args):
    # Constants
    tree_extensions = (".tree", ".tre")
    # Checks
    if args.tree.endswith(tree_extensions):
        args.base = '.'.join(args.tree.split('.')[:-1])
        if not args.output:
            args.output = args.base + '_functional_decoration.tree'            
    else:
        raise BadTreeFileException("Tree file without one of the following \
Extensions was provided: %s" % ' '.join(tree_extensions))
    return args

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''Functionally decorate a  
tree as best as possible''')
    parser.add_argument('--tree', 
                        help='Path to newick tree file',
                        required=True)
    parser.add_argument('--database', 
                        help='Path to database file - This is required \
for everything except depth clading.')
    parser.add_argument('--depth', 
                        help='Attempt to clade by depth (provide a \
percentile cutoff, or a list of percentiles (space separated)', 
                        nargs='+',
                        type = int)
    parser.add_argument('--context', 
                        help='Attempt to clade by genetic context (deprecated)', 
                        action="store_true")
    parser.add_argument('--log', 
                        help='Output logging information to this file.')
    parser.add_argument('--verbosity', 
                        help='Level of verbosity (1 - 5 - default = 4) \
5 = Very verbose, 1 = Silent',
                        type = int,
                        default = 4)
    parser.add_argument('--output', 
                        help='Output decorated tree to this file path.')
    
    args = check_args(parser.parse_args())

    if args.log:
        if os.path.isfile(args.log): raise Exception("File %s exists" % args.log)
        logging.basicConfig(filename=args.log, level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    fd = FunDec(args.base)
    
    fd.main(args.tree,
            args.output,
            args.database,
            args.depth,
            args.context)

    exit(0)
    
