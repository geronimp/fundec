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

from skbio.tree import TreeNode
from collections import Counter

# Local imports
from reroot import Rerooter

###############################################################################
############################### - Exceptions - ################################

class BadTreeFileException(Exception):
    pass

###############################################################################
################################ - Classes - ##################################


class Gff:
    def __init__(self, gff_file_obj): 
        '''
        Parse a gff file into a dictionary
        
        Parameters
        ----------
        gff_file_path: str
            Path to gff file.
        Returns
        -------
        Hash containing the parsed content of the gff file:
            {gene_id_1: annotation}
            ...
            }
        '''
        idx = 1
        self.context = {}
        self.info = {}
        for line in gff_file_obj:
            if line == "##FASTA\n": break
            elif not line.startswith('#'):
                if line.split('\t')[2] ==  "CRISPR" or \
                                line.split('\t')[2] == "direct" or \
                                line.split('\t')[2] == "repeat_region" or \
                                line.split('\t')[2] == "tandem" or \
                                line.split('\t')[2] == "" or \
                                line.split('\t')[2] == "flanking" or \
                                line.split('\t')[2] == "inverted":
                    
                    sl=line.strip().split('\t')
                    gene_start_pos=sl[3]
                    gene_finish_pos=sl[4]
                    direction = None
                    entry={"ID":None}
                                            
                else:
                    (contig, img_version, entry_type, gene_start_pos, 
                     gene_finish_pos, _, direction, _, entry) \
                                            = tuple(line.strip().split('\t')) 
                    tmp_entry = entry.split('=')[-1].translate(None, ';,-()[]')
                    entry = '='.join(entry.split('=')[:-1]+[tmp_entry])
                    
                    entry = {x.split('=')[0]:x.split('=')[1] 
                             for x in entry.replace(';;', ';').split(';')}

            
                gene_length = int(gene_finish_pos) - int(gene_start_pos)
                self.context[entry["ID"]] = idx
                self.info[idx] = [entry["ID"], 
                                  gene_length, 
                                  (True if direction == "+" else False)]
                if "product" in entry:
                    self.info[idx].append(entry['product'])
                idx+=1
                
    def surrounding(self, gene_id, distance):
        '''
        Parse out gene context of a given gene upstream and downstream by a 
        distance specified in "distance".
        
        Parameters
        ----------
        gene_id: str
            ID of gene to query
        distance: int
            The number of genes to go upstream and downstream of the
            query.
            
        Returns
        -------
        list containing genes upstream and downstream to the query by distance.
 
        '''
        logging.debug("Extracting gene context for %s" % gene_id)
        context = []
        gene_index = self.context[gene_id]+1
        for forward_idx, reverse_idx \
                        in zip(range(gene_index, (gene_index+distance+1)),
                               range((gene_index-distance-1), gene_index-1)):
            if forward_idx in self.info:
                context.append(self.info[forward_idx][0])
            if reverse_idx in self.info:
                context.append(self.info[reverse_idx][0])
        return context
    
    def identity(self, gene_id):
        gene_info = self.info[self.context[gene_id]]
        annotation = gene_info[-1]
        if type(annotation) == str:
            return annotation
        else:
            return "unclassified_protein"
        