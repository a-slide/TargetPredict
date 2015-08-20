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
from tempfile import mkdtemp
from shutil import rmtree
from itertools import islice
import shlex

# Third party package

# Local packages
from Blastn import Blastn
from MirandaHit import MirandaHit
from FileUtils import is_readable_file, is_gziped, gunzip, file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class TargetPredict (object):
    """
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "TargetPredict 0.1"
    USAGE = "Usage: %prog -s SUBJECT -q QUERY"

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

        ### Parse arguments
        opt, args = optparser.parse_args()

        # Verify the presence of options
        try:
            assert opt.subject, "Missing subject (-s) option"
            assert opt.query, "Missing query (-q) option"

        except AssertionError as E:
            print (E)
            optparser.print_help()
            sys_exit()

        ### Init a RefMasker object
        return TargetPredict (subject=opt.subject, query=opt.query)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, subject, query):
        """
        General initialization function for import and command line
        """
        print("\nInitialize TargetPredict")

        # Verify files readability
        assert is_readable_file(subject), "{} is not readable".format(subject)
        assert is_readable_file(query), "{} is not readable".format(query)

        # Create a temporary folder for working files
        self.temp_dir = mkdtemp()

        # Extract gzip in temporary files if needed
        self.subject = gunzip (subject, self.temp_dir) if is_gziped(subject) else subject
        self.query = gunzip (query, self.temp_dir) if is_gziped(query) else query

        # Define additional self variables
        self.basename = "{}_{}".format(file_basename(subject), file_basename(query))
        self.blast_hit = []
        self.miranda_hits =[]

    # Enter and exit are defined to use the context manager "with"
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        print ("Cleaning up temporary files".format(self.temp_dir))
        rmtree(self.temp_dir)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        """
        ##### BLAST prediction #####
        print ("Finding hits with BLASTN")

        # Blast using an existing wrapper/parser
        with Blastn(ref_path=self.subject) as blastn:
            self.blast_hit = blastn (query_path=self.query, task="blastn-short", evalue=1000, blastn_opt="-word_size 5")

        # Suppress hits on positive strand and sort by score
        self.blast_hit = [hit for hit in self.blast_hit if hit.strand == "-"]
        self.blast_hit.sort (key=lambda x: x.score, reverse=True)

        # Write a complete blast report
        self._write_report (self.blast_hit, "{}_raw_blast_results.csv".format(self.basename))

        ##### MIRANDA prediction + parsing #####
        print ("Finding hits with MIRANDA")

        miranda_exec = "miranda"
        miranda_score = 180
        miranda_opt = "-quiet"

        cmd = "{} {} {} {} -sc {}".format(miranda_exec, self.query, self.subject, miranda_opt, miranda_score)
        print (cmd)
        miranda_output = self.yield_cmd(cmd)

        for line in miranda_output:
            if line[0] == ">" and line[1] != ">":
                hit_split = line[1:].strip().split()
                assert len(hit_split) == 11, "Invalid miranda line: {}".format(line)
                self.miranda_hits.append (MirandaHit(*hit_split))

        self.miranda_hits.sort (key=lambda x: x.score, reverse=True)

        # Write a complete blast report
        self._write_report (self.miranda_hits, "{}_raw_miranda_results.csv".format(self.basename))

        # COMPARISON WITH GFF FILES

        # INTERSECTION (VEIN DIAGRAM)


    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def yield_cmd (self, cmd):
        """
        Decompose shell command in list of elementary elements and parse the output line by
        line using a yield statement
        """
        # Split the commands
        split_cmd = shlex.split(cmd)

        # Prepare the popen object
        proc = Popen(split_cmd, stdout=PIPE)

        # yield results line by line
        for line in proc.stdout:
            yield line


    def _write_report (self, hit_list, file_name):

        with open (file_name, "w") as fout:

            # Write table header
            fout.write ("\t".join(hit_list[0].report().keys())+"\n")

            # Write values
            for hit in hit_list:
                fout.write ("\t".join([str(val) for val in hit.report().values()])+"\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    start_time = time()
    with TargetPredict.class_init() as target_predict:
        target_predict()
    print ("\nDone in {}s".format(round(time()-start_time, 3)))
