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

# Local imports
from reroot import Rerooter

###############################################################################
############################### - Exceptions - ################################

class MisformattedTable(Exception):
    pass

###############################################################################
################################ - Classes - ##################################

class ReadDatabase:
    
    _TYPE_NUMERIC = "numeric"
    _TYPE_STRING = "string"
    _TYPE_TAXONOMY = "taxonomy"
    
    def __init__(self, database, deliminator):
        logging.info("Parsing input database")
        self.database = {}
        
        self.fieldnames = database.readline()\
                                  .strip()\
                                  .split(deliminator)
        logging.debug("%i columns found" % (len(self.fieldnames)))
        database_hash = {fieldname:{} for fieldname in self.fieldnames}
        
        for line in database:    
            columns = line.strip().split(deliminator)
            if len(columns)!=len(self.fieldnames):
                raise MisformattedTable("Number of columns in database are \
not the same for all rows!")
                
            for index, column_entry in enumerate(columns):
                if index == 0:
                    row_name = column_entry
                else:
                    fieldname = self.fieldnames[index]
                    database_hash[fieldname][row_name] = column_entry
                    
        for fieldname, entries in database_hash.iteritems():
            
            if fieldname == "taxonomy":
                type = self._TYPE_TAXONOMY
            else:    
                type = self._TYPE_NUMERIC
                for index, entry in enumerate(entries.values()):
                    if not entry.isdigit():
                        type = self._TYPE_STRING
                    
                if type == self._TYPE_NUMERIC:
                    for entry, item in entries.iteritems():
                        entries[entry] = int(item) 
                        
            logging.debug("Column %s typed as %s" % (fieldname, 
                                                     type))
            
            self.database[fieldname] = DatabaseField(fieldname,
                                                     entries,
                                                     type)

    def columns(self):
        return self.fieldnames[1:]

    def field(self, field):
        return self.database[field]
    
class DatabaseField(ReadDatabase):
    
    def __init__(self, name, entries, type):
        self.name = name
        self.entries = entries
        self.values = entries.values()
        self.type = type
        
        if type == ReadDatabase._TYPE_TAXONOMY:
            for key, tax_string in self.entries.iteritems():
                dom, phy, cla, ord, fam, gen, spe = (tax_string.split('; ')
                                                     if '; ' in tax_string
                                                     else tax_string.split(';'))
                self.entries[key] = {'domain'  : dom,
                                     'phylum'  : phy,
                                     'class'   : cla,
                                     'order'   : ord,
                                     'family'  : fam,
                                     'genus'   : gen,
                                     'species' : spe}
        
        
        
        
        
        
    