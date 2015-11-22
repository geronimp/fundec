#!/usr/bin/env python
###############################################################################
#
# fundec.py <tree> <database_file> <%id cutoff>
# I should probably iclude some reguular expression recognition to make the
# counts of each annotation more robust. I imagine that as is this
# code will only cluster the major branches. I should really set the
# %ID cutoff high to avoid over clustering. I dont want that to come
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
from gff import Gff
from tree_clust import Cluster

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
        self.HEADER_DB="/srv/home/s4293811/04_gits/fundec/database/header_db.txt"
        self.IMG_GENOME_DB = "/srv/db/img/latest/genomes/"
        self.GTDB_GENOME_DB = "/srv/projects/share/for_joel_gtdbgenomes_prokka/ace_genomes_prokka/"
    # ------------------------------- Parsing ------------------------------- #
    def _parse_db_file(self, database):
        db_dict = {}
        for line in open(database):
            if line.startswith('name'): continue
            line=line.strip().split('\t')
            db_dict[line[0]]= {'name':line[0],
                               'full_name':line[1],
                               'annotation':line[2],
                               'tax_string':line[3],
                               'pfam_ids':line[4],
                               'pfam_classifications':line[5],
                               'cog_ids':line[6],
                               'cog_classifications':line[7],
                               'tigrfam_ids':line[8],
                               'tigrfam_classifications':line[9],
                               'source':line[10],
                               'Swiss_Prot_annotated':line[11],
                               'sequence_length':line[12],
                               'domain':line[13],
                               'phylum': line[14],
                               'class': line[15],
                               'order': line[16],
                               'family': line[17],
                               'genus': line[18],
                               'species': line[19],
                               'paralog': line[20],
                               'genome_id': line[21],
                               'Alignment': line[22],
                               'classificaton':['%s_homolog' % self.db_name]}
        return db_dict

    def _open_tree(self, tree_path):
        tree_obj=TreeNode.read(open(tree_path))
        return tree_obj
    
    # ---------------------------- Manipulating ----------------------------- #
    def _gather_annotations(self, tips, db):
        'Gather annotations into a list, return relative abundances of each'
        annotations = []
        for tip in tips:
            annotations.append(db[tip.name.replace(' ', '_')]['annotation'])
        counted_annotations = Counter(annotations)
        total_annotations=float(sum(counted_annotations.values()))
        counted_annotations_relab={}
        for annotation, count in counted_annotations.iteritems():
            counted_annotations_relab[annotation]=\
                                        float(count)/total_annotations
        return counted_annotations_relab

    def  _gather_ancestors(self, node):
        previous_clade = []
        for ancestor in node.ancestors():
            if ancestor.name:
                if ancestor.name.startswith('a__'):
                    previous_clade.append(ancestor.name)
        return previous_clade

    # =============================== MAIN ================================== #
    def main(self, tree_path, database_path, output_path, percentile,
             length, context_distance, pfam, classification):
        if database_path:
            database = self._parse_db_file(database_path)
        tree     = self._open_tree(tree_path)
        c = Cluster()
        Rerooter().reroot(tree)
        
        if percentile:
            logging.info("Identifying depth clusters in tree using %ith percentile\
 cutoff" % percentile)
            tree=c.depth_first_cluster(tree, percentile)
        if length:
            logging.info("Identifying length clusters in tree")
            tree=c.length_second_cluster(tree, 
                                    {name:int(database[name]\
                                                      ['sequence_length'])
                                     for name in [x.name.replace(' ', '_') 
                                                 for x in tree.tips()]})
            
        if classification:
            logging.info("Identifying consistent classification patterns in \
tree")
            annotations = {key:item['annotation']
                                for key, item in database.iteritems()}
            if len(set(annotations.values()))==1:
                logging.info("All tips share the same annotation. \
Nothing to annotate!")
            else:
                tree=c.annotation_cluster(tree, annotations, "annotation")
                
                
        if pfam:
            logging.info("Identifying pfam annotation clusters in tree")
            pfam_annotations = {key:','.join(sorted(item\
                                        ['pfam_classifications'].split(','))) 
                                for key, item in database.iteritems()}
            if len(set(pfam_annotations.values()))==1:
                logging.info("All tips share the same pfam annotations. \
Nothing to annotate!")
            else:
                tree=c.annotation_cluster(tree, pfam_annotations, "PFAM")
        

        
        if context_distance:
            logging.info("Identifying clusters in tree with different gene \
context")   
            # Hard coding paths is a no-no!
            GFF_DB = "/srv/projects/abisko/Joel/99_phd/01_Projects/08_gpkgs_for_ko_groups/data/gff_database/"
            genome_to_gff = {}
            gene_to_context_ids = {}
            gene_to_genome_file = {}
            poor_headers = {}
            for tip in tree.tips():
                name = tip.name.replace(' ', '_')
                genome_id = database[name]['genome_id']
                ############################################################
                ############################################################
                ## This part is real hacky and a result of me pulling in  ##
                ## sequences from multiple databases and in part because  ##
                ## of my poor planning.                                   ##
                if '~' in name:
                    if "~U_" in name:
                        poor_headers[name]='_'.join(name.split('~')[0]\
                                                        .split('_')[-2:])
                        gene_to_genome_file[name] = os.path.join(self.GTDB_GENOME_DB,
                                                                "%s.faa" % genome_id) 
                    else:
                        cmd = "grep '%s' %s" % (name, self.HEADER_DB)
                        
                        poor_headers[name]=subprocess\
                                            .check_output(cmd, shell=True)\
                                            .strip().split(',')[1]  
                        gene_to_genome_file[name] = os.path.join(self.IMG_GENOME_DB,
                                                                genome_id,
                                                                "%s.genes.faa" % genome_id)    

                else:
                    gene_to_genome_file[name] = os.path.join(self.IMG_GENOME_DB,
                                        genome_id,
                                        "%s.genes.faa" % genome_id)        
                ############################################################
                ############################################################   
            
                gff_path = os.path.join(GFF_DB,
                                        database[name]['source'].upper(),
                                        "%s.gff" % genome_id)
                if genome_id not in genome_to_gff:
                    gff=Gff(open(gff_path))
                    gene_to_context_ids[name] = gff\
                                        .surrounding((poor_headers[name] if '~' in name 
                                                      else name), 
                                                      1)
                    genome_to_gff[genome_id] = gff
                else:
                    gene_to_context_ids[name] = genome_to_gff[genome_id]\
                                        .surrounding((poor_headers[name] if '~' in name 
                                                      else name), 
                                                      1)
            tree=c.synteny_cluster(tree, gene_to_context_ids, 
                                   gene_to_genome_file, self.db_name,
                                   genome_to_gff, database)
            
        # Old code. Not sure if redundant.
        #for node in tree.preorder():
        #    tip_annotations = self._gather_annotations(node.tips(),
        #                                               database)
        #    if any([x for x in tip_annotations.values() if x > cutoff]):
        #        previous_clades = self._gather_ancestors(node)
        #        key, _ = max(tip_annotations.iteritems(), key=lambda x:x[1])
        #        new_classification = "a__%s" % key
        #        if new_classification not in previous_clades:
        #            node.name = new_classification
        logging.info("Finished clustering! Writing tree to file: %s" \
                                                            % (output_path))
        tree.write(output_path, format = "newick")

###############################################################################
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
    if not args.database:
        if not args.length and not args.pfam and not args.context \
                           and not args.annotation:
            pass
        else:
            logging.error("A database must be provided to use any of the \
annotation options except depth.")
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
percentile cutoff - recommended = 10)', 
                        type = int)
    parser.add_argument('--length', 
                        help='Attempt to clade by length', 
                        action="store_true")
    parser.add_argument('--pfam', 
                        help='Attempt to clade by pfam annotation', 
                        action="store_true")
    parser.add_argument('--annotation', 
                        help='Attempt to clade by classification', 
                        action="store_true")
    parser.add_argument('--context', 
                        help='Attempt to clade by genetic context', 
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
            args.database,
            args.output,
            args.depth,
            args.length,
            args.context,
            args.pfam,
            args.annotation)

    exit(0)
    
