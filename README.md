# fundec
Tree decoration - made fun!
This is a small tool for decorating newick format phylogenetic trees with customisable information

Fundec requires a rooted phylogenetic tree in newick format as input. You can use several approaches to decoration:
* Depth clustering using the partition method modified from:
```
        Prosperi, M.C.F., et al. A novel methodology for large-scale phylogeny
        partition. Nat. Commun. 2:321 doi: 10.1038/ncomms1325 (2011).
        http://www.nature.com/ncomms/journal/v2/n5/full/ncomms1325.html
```
* Taxonomy decoration
* Decorate with information from an input table (Can be numeric or string based)
* Any combination of the above

### Usage
Any of the following grouping methods can be used in combination with one another:
####Depth
To cluster using tree depth, use the --depth flag, followed by the percentile cutoffs you wish to cluster at. The higher the number, the larger the clades. You can specify more than one percentile:
```
fundec.py --depth 5 10 20 --tree my.tree --output my_decorated.tree
```
####Decoration with database file
A table can be given to fundec with the sequence name as the first column, followed by columns filled with information you wish to decorate with:

There are three types of information you can decorate with:
*''numeric'' numeric based decoration will occur if all rows in a column are numbers. Approach will iterate through the tree is post order, and use a basic student's t-test to identify clades with significantly different values to that of its closest sibling clade.
*''string'' String based decoration will iterate through the tree and cluster node that have a > 80% consensus of all tips from that node.
*''taxonomy'' Taxonomy based decoration is a special instance of string based decoration. A column with the header 'taxonomy' filled with taxonomy strings will automatically decorate the tree with taxonomy. The approach used by fundec is simplistic in that it will not cluster clades with any inconsistent taxonomy. This means it isfar more strict than software such as [Tax2Tree](https://github.com/biocore/tax2tree). Expect clades to be broken.

Provide the table as either a .csv or a .tsv file to the --database flag as in the following example:

```
fundec.py --database my_database.tsv --tree my.tree --output my_output.tree
```


### Contact
If you have any further comments, complaints or recommendations about fundec, please raise an issue on this page.
Software by [Joel A. Boyd](http://ecogenomic.org/users/joel-boyd) (@geronimp) at the [Australian Centre for Ecogenomics](http://ecogenomic.org)
Released under GPL3 - see LICENCE.txt for licensing details
