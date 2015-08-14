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
from sys import exit as sys_exit
from time import time
from subprocess import Popen, PIPE
#from gzip import open as gopen
from collections import OrderedDict
from datetime import datetime

# Third party package

# Local packages
from Blastn import Blastn


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

        # Verify the presence of options
        try:
            assert opt.subject, "Missing subject (-s) option"
            assert opt.gff, "Missing query (-q) option"

        except AssertionError as E:
            print (E)
            optparser.print_help()
            sys_exit()


        ### Init a RefMasker object
        return TargetPredict (subject=opt.subject, query=opt.query) #, gff=opt.gff)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, subject, query):#, gff=None):
        """
        General initialization function for import and command line
        """
        print("\nInitialize TargetPredict")

         ### Storing Variables
        self.subject = subject
        self.query = query
        #self.gff = gff
        self.basename = "{}_{}".format(self.file_basename(subject), self.file_basename(query))
        self.blast_hit = []
        self.miranda_hits =[]


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        """
        start_time = time()

        # BLAST prediction
        print ("Finding hits with blast")

        # Using an existing wrapper/parser
        with Blastn(ref_path=self.subject) as blastn:
            self.blast_hit = blastn( query_path=self.query, task="blastn-short", evalue=1000, blastn_opt="-word_size 5")

            # filter hit on negative strand only
            # Select the n better alignements

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
