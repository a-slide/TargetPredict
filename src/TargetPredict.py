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
import optparse
from sys import exit as sys_exit
from time import time
from subprocess import Popen, PIPE
from collections import OrderedDict
from datetime import datetime
from tempfile import mkdtemp
from shutil import rmtree
import shlex
from multiprocessing import cpu_count

# Local packages
from BlastHit import BlastHit
from MirandaHit import MirandaHit
from FileUtils import is_readable_file, is_gziped, gunzip, file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class TargetPredict (object):
    """
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "TargetPredict 0.1"
    USAGE = "Usage: %prog -s SUBJECT -q QUERY [...]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-s', '--subject', dest="subject",
            help= "Path to the subject in fasta")
        optparser.add_option('-q', '--query', dest="query",
            help= "Path to the query in fasta")
        default = "-task blastn-short -evalue 100000000 -word_size 5 -max_target_seqs 10000"
        optparser.add_option('--blastn_opt', dest="blastn_opt", default=default,
            help= "Options to run blastn (default='{}')".format(default))
        default = "-sc 160"
        optparser.add_option('--miranda_opt', dest="miranda_opt", default=default,
            help= "Options to run miranda (default='{}')".format(default))

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
        return TargetPredict (
            subject = opt.subject,
            query = opt.query,
            blastn_opt = opt.blastn_opt,
            miranda_opt = opt.miranda_opt)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, subject, query, blastn_opt="", miranda_opt=""):
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
        self.original_subject = subject
        self.orignal_query = query
        self.subject = gunzip (subject, self.temp_dir) if is_gziped(subject) else subject
        self.query = gunzip (query, self.temp_dir) if is_gziped(query) else query
        self.blastn_opt = blastn_opt
        self.miranda_opt = miranda_opt

        # Define additional self variables
        self.basename = "{}_{}".format(file_basename(subject), file_basename(query))

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

        print ("\nFinding hits with BLASTN")

        blast_hits =[]

        # Create blastn database
        print ("  Create a blast database")
        db_path = "{}/{}".format(self.temp_dir, file_basename(self.subject))
        makeblastdb_cmd = "makeblastdb -dbtype nucl -input_type fasta -in {} -out {}".format(self.subject, db_path)
        print ("\t"+makeblastdb_cmd)
        makeblastdb_out = self._yield_cmd(makeblastdb_cmd)

        with open (file_basename(self.subject)+"_makeblastdb.log", "w") as fout:
            for line in makeblastdb_out:
                fout.write(line)

        # Perform blast
        print ("  Run blastn")
        blastn_cmd = "blastn {} -num_threads {} -outfmt \"6 std qseq\" -dust no -query {} -db {}".format(self.blastn_opt, cpu_count(), self.query, db_path)
        print ("\t"+blastn_cmd)

        blastn_out = self._yield_cmd(blastn_cmd)
        for line in blastn_out:
            hit_split = line.strip().split()
            assert len(hit_split) == 13, "Invalid blast line: {}".format(line)
            blast_hits.append(BlastHit(*hit_split))

        # Suppress hits on positive strand, keep best hit per subject only and sort by score
        print ("  Process hits")
        blast_hits = [hit for hit in blast_hits if hit.strand == "-"]

        blast_hit_dict = {}
        for hit in blast_hits:
            if hit.s_id in blast_hit_dict:
                if hit.score > blast_hit_dict[hit.s_id].score:
                    blast_hit_dict[hit.s_id]=hit
            else:
                blast_hit_dict[hit.s_id]=hit
        blast_hits = blast_hit_dict.values()

        blast_hits.sort (key=lambda x: x.score, reverse=True)

        # Write a complete blast report
        print ("  Write a blast report")
        self._write_report (blast_hits, "{}_raw_blast_results.csv".format(self.basename))

        ##### MIRANDA prediction#####

        print ("\nFinding hits with MIRANDA")

        miranda_hits =[]

        # Run miranda and parse output
        print ("  Run miranda")
        miranda_cmd = "miranda {} {} -quiet {}".format(self.query, self.subject, self.miranda_opt)
        print ("\t"+miranda_cmd)

        miranda_output = self._yield_cmd(miranda_cmd)

        for line in miranda_output:
            if line[0] == ">" and line[1] != ">":
                hit_split = line[1:].strip().split()
                assert len(hit_split) == 11, "Invalid miranda line: {}".format(line)
                miranda_hits.append (MirandaHit(*hit_split))

        # Keep best hit per subject only and sort hits by score
        print ("  Process hits")
        miranda_hit_dict = {}
        for hit in miranda_hits:
            if hit.s_id in miranda_hit_dict:
                if hit.score > miranda_hit_dict[hit.s_id].score:
                    miranda_hit_dict[hit.s_id]=hit
            else:
                miranda_hit_dict[hit.s_id]=hit
        miranda_hits = miranda_hit_dict.values()

        miranda_hits.sort (key=lambda x: x.score, reverse=True)

        # Write a complete miranda report
        print ("  Write a Miranda report")
        self._write_report (miranda_hits, "{}_raw_miranda_results.csv".format(self.basename))

        # Write a report
        report_out = self.basename+".report.txt"
        print ("\nGenerate a summary report")
        with open (report_out, "w") as fout:
            fout.write ("Program {}\tDate {}\n".format(self.VERSION,str(datetime.today())))
            fout.write ("\n### OPTIONS ###\n")
            fout.write ("Subject fasta file\t{}\n".format(self.original_subject))
            fout.write ("Query fasta file\t{}\n".format(self.orignal_query))
            fout.write ("Makeblastdb command\t{}\n".format(makeblastdb_cmd))
            fout.write ("Blastn command\t{}\n".format(blastn_cmd))
            fout.write ("Miranda command\t{}\n".format(miranda_cmd))
            fout.write ("\n### COUNTS ###\n")
            fout.write ("Blast Hits found\t{}\n".format(len(blast_hits)))
            fout.write ("Miranda Hits found\t{}\n".format(len(miranda_hits)))

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _yield_cmd (self, cmd):
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
