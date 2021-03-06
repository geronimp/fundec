#!/usr/bin/env python
###############################################################################
#
# fundec.py <tree> <database_file> <%id cutoff>
# I should probably include some regular expression recognition to make the
# counts of each annotation more robust. I imagine that as is this
# code will only cluster the major branches. I should really set the
# %ID cutoff high to avoid over clustering. I dont want that to come
# back to bite.
#
# Come to think of it, would some kind of modification of Levenshtein 
# distance work here?
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
import logging
import tempfile
import subprocess
import numpy as np
import os
import re
import random

from scipy import stats
from itertools import combinations
from collections import Counter
from skbio.tree import TreeNode

###############################################################################
############################### - Exceptions - ################################

class MalformedTreeException(Exception):
    pass

###############################################################################
################################ - Classes - ##################################

class Cluster:

    def _node_dist(self, node):
        '''
        Returns a list of distances tip-to-tip from the node provided. 
        
        Parameters
        ----------
        node: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
        Returns
        -------
        Array of distances from tip-to-tip
        '''
        logging.debug("Calculating tip-to-tip distances for node: %s" \
                                                                % node.name)
        distances=[]
        if len(list(node.tips())) > 200:
            
            for tip_a, tip_b in list(random.sample(list(combinations(node.tips(), 2)), 6000)):
                distances.append(tip_a.distance(tip_b))
            distances=np.array(sorted(distances))
        else:
            for tip_a, tip_b in combinations(node.tips(), 2):
                distances.append(tip_a.distance(tip_b))
            distances=np.array(sorted(distances))
        return distances
    
    def _rename(self, node, name):
        
        if node.name:
            try: 
                float(node.name)
                node.name = "%s:%s" % (node.name,
                                       name)
            except:
                node.name = "%s; %s" % (node.name,
                                        name)
        else:
            node.name = name
                   
    def depth_first_cluster(self, tree, percentile, index):
        '''
        Attempt to cluster tree with nodes of tip-to-tip distrubution <
        an nth percentile cutoff of the whole-tree distance distribution. 
        A better description can be found in the citation below.
        
        Parameters
        ----------
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
                
        percentile: float
            The percentile cutoff to use to determine the cutoff from clading
            from a given node.
        
        index: int
            Number indicating the current iteration round of clustering. This 
            number assists the distinction of depth clusters.
        
        Clustering method modified from Prosperi et al method:
        Prosperi, M.C.F., et al. A novel methodology for large-scale phylogeny 
        partition. Nat. Commun. 2:321 doi: 10.1038/ncomms1325 (2011).
        
        http://www.nature.com/ncomms/journal/v2/n5/full/ncomms1325.html
        '''
            
        cluster_count = 1
        clustered = []
        
        if tree.is_root():
            logging.debug("Calculating %ith percentile cutoff from root" \
                                                    % (percentile))
            whole_tree_distribution = self._node_dist(tree)
            
            cutoff = np.percentile(whole_tree_distribution, percentile)
            logging.debug("Cutoff (%ith percentile): %f" % (percentile,
                                                            cutoff))
        else:
            raise MalformedTreeException('''Could not determine the \ 
distance distribution of whole tree. Did you provide the whole tree or just a \
node?''')
        
        for node in tree.preorder():
            if node in clustered:
                continue
            elif node.is_tip():
                continue
            else:
                node_distribution = self._node_dist(node)
                median=np.median(node_distribution)
                logging.debug("Median of node: %f" % median)
                if median <= cutoff:
                    logging.debug("Cluster found!")
                    self._rename(node, "t__%i-%i" % (index, cluster_count))
                    cluster_count+=1
                    for descenent in node.traverse():
                        clustered.append(descenent)
        logging.info("%i depth cluster(s) found in tree" % (cluster_count-1))
        return tree
    
    def numeric_decoration(self, tree, leaf_lengths, prefix):
        '''
        Attempt to cluster tree based on gene length within clades
        
        Parameters
        ----------
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
        leaf_lengths: hash
            dictionary with each leaf name as the key, and its respective gene
            length (int) the entry. As follows:
                {
                gene_1: 456
                gene_2: 483
                ...
                }
        prefix: str
            string to distinguish different annotations when decorating
        '''
        
        cluster_count = 1
        clustered = []
        
        for node in tree.postorder():
            if node.is_root():
                continue
            elif any([x for x in node.traverse() if x in clustered]):
                continue
            elif node in clustered:
                continue
            elif node.is_tip():
                continue
            else:
                tip_names = set([tip.name.replace(' ', '_') 
                                 for tip in node.tips()])
                
                in_node = [item for key, item in leaf_lengths.iteritems()
                           if key in tip_names]
                if len(node.siblings()) > 1:
                    raise MalformedTreeException("Node encountered with more \
than one sibling!")
                    exit(1)
                else:
                    sibling = node.siblings()[0]
                    if sibling.is_tip():
                        continue    
                    else:
                        out_node = [leaf_lengths[x.name\
                                                  .replace(' ', '_')] 
                                    for x in sibling.tips()]
                    t, p = stats.ttest_ind(out_node, in_node)     
                    in_mean = int(sum(in_node)/len(in_node))           
                    out_mean = int(sum(out_node)/len(out_node))
                    if p < 0.0001:
                        self._rename(node, "%s__%i" % (prefix, in_mean))
                        cluster_count += 1
                        for descenent in node.traverse():
                            clustered.append(descenent)
                        clustered.append(node)
            
                    else:
                        cluster_count += 1
                        pre_clustered=False
                        for x in node.parent.non_tips():
                            if x.name:
                                if ("%s__" % prefix) in x.name:
                                    pre_clustered=True
                        if pre_clustered:
                            mean = int(sum(in_node)/len(in_node))
                            self._rename(node, "%s__%i" % (prefix, mean))
                            
                            for descenent in node.non_tips():
                                clustered.append(descenent)
                            clustered.append(node)  
                        else:
                            mean = int(sum(out_node + in_node)/len(out_node + in_node))
                            self._rename(node.parent, "%s__%i" % (prefix, mean))
                            
                            for descenent in node.parent.non_tips():
                                clustered.append(descenent)
                            clustered.append(node)    
                             
                            
        logging.info("%i length cluster(s) found in tree" % (cluster_count-1))
        return tree
    
    def string_decoration(self, tree, annotations, prefix):
        '''
        Attempt to cluster tree based on pfam annotations awarded to eahc leaf
        
        Parameters
        ----------
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
        annotations: hash
            dictionary with each leaf name as the key, and its respective 
            annotation as the entry. As follows:
                {
                gene_1: 'MCR_alpha,MCR_alpha_N'
                gene_2: 'MCR_alpha'
                ...
                }
        prefix: str
            string to distinguish different annotations when decorating
        '''
        cluster_count = 1
        clustered = []
        for node in tree.preorder():
            if node.is_root():
                continue
            elif node in clustered:
                continue
            elif node.is_tip():
                continue
            else:
                tip_names = set([tip.name.replace(' ', '_') 
                                 for tip in node.tips()])

                in_node_pfam = [annotations[x] for x in tip_names]

                total=float(len(in_node_pfam))
                
                in_node_pfam_fraction = {float(item)/total:key 
                                         for key, item in Counter(in_node_pfam)\
                                                                .iteritems()}

                most_common = max(in_node_pfam_fraction.keys())
                
                if most_common > 0.90:
                    logging.debug("Cluster found!")
                    self._rename(node, "%s__%s" \
                                        % (prefix, 
                                           in_node_pfam_fraction[most_common]))
                    cluster_count += 1
                    for descenent in node.traverse():
                        clustered.append(descenent)
        logging.info("%i %s cluster(s) found in tree" % (cluster_count-1, 
                                                         prefix))
        return tree                    
                    
    def taxonomy_decoration(self, tree, gene_to_taxonomy):
        '''
        Partition the tree into strict taxonomy clades. This code currently 
        does not allow for inconsistency within a clade of any fashion, but 
        should be flexible enough to integrate easily should the need arise.
        As a result, taxonomy decoration can be quite messy with this code,
        particularly if you have many HGT events. Because this causes 
        paraphyletic groups in what is supposed to be monophyletic, a simple 
        tally system distinguished taxonomic groupings.
        
        Parameters
        ----------
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
        
        gene_to_taxonomy: dict
            In the following format:
                
                {tip_name_1: 'class': 'c__Gammaproteobacteria',
                             'domain': 'd__Bacteria',
                             'family': 'f__Pasteurellaceae',
                             'genus': 'g__Actinobacillus_Mannheimia',
                             'order': 'o__Pasteurellales',
                             'phylum': 'p__Proteobacteria',
                             'species': 's__varigena'}
                 .... etc
                 }
                 
        Returns
        -------
        skbio TreeNode object'''
        
        taxonomy_array = ["domain", "phylum", "class", "order", "family", 
                          "genus", "species"]
        assigned_nodes = set()
        
        for node in tree.preorder():
            if node not in assigned_nodes:
                
                tax_string_array = []
                
                for rank in taxonomy_array:
                    rank_tax = \
                        Counter([gene_to_taxonomy[tip\
                                                    .name\
                                                    .replace(' ','_')][rank] 
                                 for tip in node.tips()])
                        
                    consistent = (True if len(rank_tax) == 1 else False)
                    
                    if consistent:
                        tax_string_array.append(rank_tax.keys()[0])
                    else:
                        continue
            tax_string_array = [x for x in tax_string_array if len(x)>3]
            if any(tax_string_array):
                
                index = 0
                for anc in node.ancestors():
                    try:
                        index+=anc.tax
                        if "d__" in anc.name:break
                    except:
                        continue
                tax_string_array = tax_string_array[index:]
                
                if len(tax_string_array) >0:
                    self._rename(node, '; '.join(tax_string_array))
                    node.tax = len(tax_string_array)
            
        return tree
                    
                    
                    
                    
                    