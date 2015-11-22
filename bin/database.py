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

name
full_name
Classification
tax_string
pfam_ids
pfam_classifications
cog_ids
cog_classifications
tigrfam_ids
tigrfam_classifications
source
Swiss_Prot_annotated
sequence_length
domain
phylum
class
order
family
genus
species
paralog
genome_id
Alignment

class TreeTips:
    
    def __init__(self, db_path):
        db_dict = {}
        for line in open(db_path):
            if line.startswith('name'): continue
            line=line.strip().split('\t')
            db_dict[line[0]]= {'name':line[0],
                               'full_name':line[1],
                               'Classification':line[2],
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
                               'phylum':line[14],
                               'class':line[15],
                               'order':line[16],
                               'family':line[17],
                               'genus':line[18],
                               'species':line[19],
                               'paralog':line[20],
                               'genome_id':line[21],
                               'Alignment':line[22]}
            
        return db_dict
    # ---------------------------- Manipulating ----------------------------- #


    # =============================== MAIN ================================== #



