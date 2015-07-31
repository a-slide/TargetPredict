#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    TargetPredict
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library packages
#from os import path, remove
import optparse
import sys
from time import time
from subprocess import Popen, PIPE
#from gzip import open as gopen
from collections import OrderedDict
from datetime import datetime

# Third party package

# Local packages
from pyBlast.Blastn import Blastn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class TargetPredict (object):
    """
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "TargetPredict 0.1"
    USAGE = "Usage: %prog -s SUBJECT -q QUERY -a ANNOTATION -o OFFSET"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        print("\nParse commande line arguments")

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-s', '--subject', dest="subject",
            help= "Path to the subject in fasta (ungziped)")
        optparser.add_option('-q', '--query', dest="query",
            help= "Path to the query in fasta (ungziped)")
        optparser.add_option('-g', '--gff', dest="gff",
            help= "Path to a gff file containing exon annotation")
        optparser.add_option('-o', '--offset', dest="offset", default = 0,
            help= "Number of bases before and after exons were to look for query hit")

        ### Parse arguments
        opt, args = optparser.parse_args()

        ### Init a RefMasker object
        return TargetPredict (subject=opt.subject, query=opt.query, gff=opt.gff, offset=opt.offset)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, subject=None, query=None, gff=None, offset=0):
        """
        General initialization function for import and command line
        """
        print("\nInitialize TargetPredict")

        ### Verify
        assert subject, "A path to a subject file is mandatory"
        assert query, "A path to a query file is mandatory"
        assert gff, "A path to a gff file is mandatory"

        ### Storing Variables
        self.subject = subject
        self.query = query
        self.gff = gff
        self.offset = int(offset)
        self.basename = "{}_{}".format(self.file_basename(subject), self.file_basename(query))


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        """
        start_time = time()
        print self.__dict__

        # BLAST prediction
        print ("Finding hits with blast")

        # Using an existing wrapper/parser
        with Blastn(ref_path=self.subject) as blastn:
            hit_list = blastn(
                query_path=self.query,
                task="blastn-short",
                evalue=20,
                blastn_opt="-word_size 5")

        print ("Found {} hits".format(len(hit_list)))

        # BLAT prediction + parsing

        #cmd = "{} > {}".format ()
        #proc = Popen(cmd, shell=True)
        #proc.communicate()[0]

        # MIRANDA prediction + parsing

        #cmd = "{} > {}".format ()
        #proc = Popen(cmd, shell=True)
        #proc.communicate()[0]

        # COMPARISON WITH GFF FILES

        # INTERSECTION (VEIN DIAGRAM)

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))
        return(0)

    def file_basename (self, path):
        """ return The basename of a file without folder location and extension """
        return path.rpartition('/')[2].partition('.')[0]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    target_predict = TargetPredict.class_init()
    target_predict()
