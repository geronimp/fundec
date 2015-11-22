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
    
    def _self_blast(self, reads):
        hit_idx_hash = {}
        
        logging.debug("Creating blast database") 
        cmd = 'makeblastdb -in %s -dbtype prot' % reads
        subprocess.call(cmd, shell=True)
        logging.debug("Self blasting")
        cmd = "blastp -db %s -query %s -out /dev/stdout -outfmt \"6 qseqid sseqid bitscore\" " \
                                                % (reads, reads)
        idx_cluster_hash = {}            
        idx = 0

        logging.debug("Parsing results (single linkage clustering)")
        for res in subprocess.check_output(cmd, shell=True).strip().split('\n'):
            query, hit, bit = res.split('\t')
            if query==hit:continue
            if float(bit)>100:
                if query in hit_idx_hash:
                    if hit in hit_idx_hash:
                        if hit_idx_hash[query] == hit_idx_hash[hit]:
                            pass
                        else:
                            try:
                                l1=idx_cluster_hash[hit_idx_hash[query]]
                                l2=idx_cluster_hash[hit_idx_hash[hit]]
                                new_l = list(set((l1 + l2)))
                                idx_cluster_hash[idx]=new_l
                                
                                del idx_cluster_hash[hit_idx_hash[query]]
                                del idx_cluster_hash[hit_idx_hash[hit]]
                                for i in new_l:
                                    hit_idx_hash[i]=idx
                                idx+=1
                            except:
                                import IPython ; IPython.embed()
                    else:
                        idx_cluster_hash[hit_idx_hash[query]].append(hit)
                        hit_idx_hash.update({hit:hit_idx_hash[query]})
                else:
                    if hit in hit_idx_hash:
                        idx_cluster_hash[hit_idx_hash[hit]].append(query)
                        hit_idx_hash.update({query:hit_idx_hash[hit]})
                    else:
                        idx_cluster_hash[idx] = [query, hit]
                        hit_idx_hash.update({query:idx, hit:idx})
                        idx+=1  
               
        return idx_cluster_hash, hit_idx_hash
            
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
        
    def depth_first_cluster(self, tree, percentile):
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
        
        Clustering method modified from Prosperi et al method:
        Prosperi, M.C.F., et al. A novel methodology for large-scale phylogeny 
        partition. Nat. Commun. 2:321 doi: 10.1038/ncomms1325 (2011).
        
        http://www.nature.com/ncomms/journal/v2/n5/full/ncomms1325.html
        '''
        
        tree = self._strip_tree(tree)
        
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
                    node.name = "depth_cluster_%i" % cluster_count
                    cluster_count+=1
                    for descenent in node.traverse():
                        clustered.append(descenent)
        logging.info("%i depth cluster(s) found in tree" % (cluster_count-1))
        return tree
    
    def length_second_cluster(self, tree, leaf_lengths):
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
                
                in_node = [item for key, item in leaf_lengths.iteritems()
                           if key in tip_names]
                out_node = [item for key, item in leaf_lengths.iteritems()
                           if key not in tip_names]
                t, p = stats.ttest_ind(out_node, in_node)     
                mean = int(sum(in_node)/len(in_node))           
                if p < 0.001:
                    logging.debug("Cluster found!")
                    if node.name:
                        node.name = node.name +  "_length_%i" % mean
                    else:
                        node.name = "length_cluster_%i" % mean
                    cluster_count += 1
                    for descenent in node.traverse():
                        clustered.append(descenent)
        logging.info("%i length cluster(s) found in tree" % (cluster_count-1))
        return tree
    
    def annotation_cluster(self, tree, annotations, suffix):
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
        suffix: str
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
                
                if most_common > 0.60:
                    logging.debug("Cluster found!")
                    if node.name:
                        node.name = node.name +  "_%s_%s" \
                                        % (suffix,
                                           in_node_pfam_fraction[most_common])
                    else:
                        node.name = "%s_cluster_%s" \
                                        % (suffix,
                                           in_node_pfam_fraction[most_common])
                    cluster_count += 1
                    for descenent in node.traverse():
                        clustered.append(descenent)
        logging.info("%i %s cluster(s) found in tree" % (cluster_count-1, 
                                                         suffix))
        return tree
    
    def synteny_cluster(self, tree, gene_to_context, gene_to_genome_file, 
                        db_name, genome_to_gff, database):
        '''
        Attempt to cluster tree by degree of synteny
        '''
        context_to_gene = {}
        gene_to_annotations={}
        context_reads = db_name+'_context.faa'

        for gene, group in gene_to_context.iteritems():
            group = [x for x in group if x]
            if any(group):
                logging.debug("Extracting genes surrounding %s" % gene)
                cmd = "fxtract -H -X -f /dev/stdin %s >> %s" \
                                                % (gene_to_genome_file[gene],
                                                   context_reads)
                
                process = subprocess.Popen(["bash", "-c", cmd], 
                                           stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE)
    
                process.communicate('\n'.join(group))
                for read in group:
                    context_to_gene[read]=gene
                    
            
            
        blast_clusters, hit_idx_hash=self._self_blast(context_reads)
        clusters = blast_clusters
        os.remove(context_reads)
        
        logging.debug("Adding singletons")
        for surrounding_group in [x for x in gene_to_context.values()]:
            for read_id in surrounding_group:
                if read_id not in hit_idx_hash:
                    idx = 1+max(clusters)
                    clusters[idx]=[read_id]

                
        logging.debug("Renaming clusters according to most abundant annotation")
        ids = list(clusters.keys())
        for key in ids:
            group = clusters[key]
            annotations = []
            for gene in group:
                if gene:
                    gff=genome_to_gff[database[context_to_gene[gene]]['genome_id']]
                    annotations.append(gff.identity(gene))
            try:
                total = len(group)
                fraction_annotation = {float(item)/total:key 
                         for key, item in Counter(annotations)\
                                                .iteritems()}
            except:
                import IPython ; IPython.embed()
            annotation = fraction_annotation[max(fraction_annotation)]
            
            clusters[annotation] = clusters.pop(key)

        logging.debug("Assigning tips with synteny information")
        for annotation, gene_group in clusters.iteritems():
            gene_group = [x for x in gene_group if x]
            
            if any(gene_group):
                
                for gene in gene_group:
                    if context_to_gene[gene] not in gene_to_annotations:
                        gene_to_annotations[context_to_gene[gene]] = [annotation]  
                    else:
                        gene_to_annotations[context_to_gene[gene]] = \
                                        sorted(gene_to_annotations[context_to_gene[gene]] + [annotation])
            else:
                gene_to_annotations[context_to_gene[gene]]=["no_annotated_genes"]
        for key, entry in gene_to_annotations.iteritems():
            gene_to_annotations[key] = ','.join(entry)

        return self.annotation_cluster(tree, gene_to_annotations, "synteny")
                    
                    
                    
                    
                    
                    
                    
                    