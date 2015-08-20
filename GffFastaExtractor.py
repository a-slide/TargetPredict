#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library packages
from gzip import open as gopen
from urllib import unquote
from collections import OrderedDict
from sys import exit as sys_exit
import optparse
from time import time
from datetime import datetime

# Third party packages
from Bio import SeqIO

# Local packages
from GffParser import GffParser

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GffFastaExtractor (object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "GffFastaExtractor 0.1"
    USAGE = "{}\n{}\n{}\n{}\n".format(
        "Usage: %prog -f FASTA_file -g GFF_file [...]",
        "Output a fasta file with sequences corresponding to the regions indicated in the GFF_file",
        "Features in the gff files can be restricted to a list of chromosomes and a list of types",
        "Features can be extended in 5' and 3' and overlapping features can be merged" )

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-f', '--fasta', dest="fasta",
            help= "Path to the fasta file contaning the complete genome (can be gzipped)")
        optparser.add_option('-g', '--gff', dest="gff",
            help= "Path to the gff file containing annotations of the genome file (can be gzipped)")
        optparser.add_option('-o', '--offset', dest="offset", default = 0,
            help= "Bases to extract before and after the feature (default 0)")
        optparser.add_option('--fusion', dest="fusion", action='store_true', default = False,
            help= "Fuse overlapping features in a meta-feature (default False)")
        optparser.add_option('--output_gff', dest="output_gff", action='store_true', default = False,
            help= "Output the gff file corresponding to the extracted features sequences (default False)")
        optparser.add_option('--features', dest="features", default = "exon",
            help= "Restrict extraction to a list of features. The list must be SPACE separated and quoted (default exon)")
        optparser.add_option('--chromosomes', dest="chromosomes", default = "",
            help= "Restrict extraction to a list of chromosomes. The list must be SPACE separated and quoted (default all)")

        ### Parse arguments
        opt, args = optparser.parse_args()

        # Verify the presence of mandatory options
        try:
            assert opt.fasta, "Missing fasta (-f) option"
            assert opt.gff, "Missing gff (-g) option"

        except AssertionError as E:
            print (E)
            optparser.print_help()
            sys_exit()

        ### Init a RefMasker object
        return GffFastaExtractor (
            fasta = opt.fasta,
            gff = opt.gff,
            offset = opt.offset,
            fusion = opt.fusion,
            output_gff = opt.output_gff,
            features = opt.features.split(),
            chromosomes = opt.chromosomes.split())

    ##### FONDAMENTAL METHODS #####

    def __init__ (self, fasta, gff, offset=0, fusion=False, output_gff=False, features=[], chromosomes=[]):
        """Init fonction parsing fasta and gff files"""

        print("Initialize GffFastaExtractor")

        self.fasta = fasta
        self.gff = gff
        self.offset = int(offset)
        self.fusion = fusion
        self.output_gff = output_gff
        self.features = features
        self.chromosomes = chromosomes
        self.separator = "|" # substitution separator

        self.out_name = "{}_{}_Offset-{}_Features-{}_Chr-{}_{}".format (
            self.fasta.rpartition('/')[2].partition('.')[0],
            self.gff.rpartition('/')[2].partition('.')[0],
            self.offset,
            "-".join(self.features) if self.features else "all",
            "-".join(self.chromosomes) if self.chromosomes else "all",
            "fused" if self.fusion else "not_fused")

        # Parse the fasta sequence and store in a dict
        print("\n  Parsing fasta file")
        openFunc = gopen if self.fasta.endswith(".gz") else open
        with openFunc (self.fasta, "r") as fasta_in:
            self.seq_dict = OrderedDict()
            all_seq = valid_seq = 0
            for seq in SeqIO.parse(fasta_in, "fasta"):
                all_seq += 1
                if not chromosomes or seq.id in chromosomes:
                    valid_seq += 1
                    self.seq_dict[seq.id]= seq
        print ("    Extract {} sequence(s) out of {}".format(
            valid_seq,
            all_seq))

        # Parse the gff file with GffParser
        print("\n  Parsing gff file")
        self.gff_parser = GffParser (
            gff_file=gff,
            offset=self.offset,
            fusion=self.fusion,
            features=self.features,
            chromosomes=self.chromosomes)

        print ("    Extract {} valid feature(s) out of {} from {} valid sequence(s) out of {}".format (
            self.gff_parser.valid_features,
            self.gff_parser.all_features,
            self.gff_parser.valid_seq,
            self.gff_parser.all_seq ))
        if fusion:
            print ("    {} out of {} features remaining after fusion".format (
                self.gff_parser.fused_features,
                self.gff_parser.valid_features))

    ##### PUBLIC METHODS #####

    def __call__ (self):
        """Launch the extraction of features """

        print("\nExtract feature sequences")

        # Write the fasta file containing the sequences of the selected features
        fasta_out = self.out_name+".fa.gz"
        print("\n  Write fasta output")
        with gopen (fasta_out, "w") as fout:
            for seq_id, gff_sequence in self.gff_parser.gff_dict.items():

                assert seq_id in self.seq_dict, "fasta and gff are incompatible: {} not found in fasta".format(seq_id)

                print ("    Extracting features from sequence {}".format(seq_id))
                for feature in gff_sequence.features:

                    fout.write(">{}\n{}\n".format(
                            str(feature).replace("\t", self.separator).replace(" ", "_"),
                            self.extract_seq(seq_id, feature.start, feature.end, feature.strand)))

        # Write the gff file containing the selected features if required
        if self.output_gff:
            gff_out = self.out_name+".gff.gz"
            print("\n  Write gff output")
            with gopen (self.out_name+".gff.gz", "w") as fout:
                fout.write(str(self.gff_parser))

        # Write a report
        report_out = self.out_name+".report.txt"
        print ("\n  Generate a summary report")
        with open (report_out, "w") as fout:
            fout.write ("Program {}\tDate {}\n".format(self.VERSION,str(datetime.today())))
            fout.write ("\n### OPTIONS ###\n")
            fout.write ("Original fasta\t{}\n".format(self.fasta))
            fout.write ("Original gff\t{}\n".format(self.gff))
            fout.write ("Offset\t{}\n".format(self.offset))
            fout.write ("Fusion\t{}\n".format(self.fusion))
            fout.write ("Output gff\t{}\n".format(self.output_gff))
            fout.write ("Restricted features\t{}\n".format("\t".join(self.features)))
            fout.write ("Restricted chromosomes\t{}\n".format("\t".join(self.chromosomes)))
            fout.write ("\n### COUNTS ###\n")
            fout.write ("Sequence(s) in gff file\t{}\n".format(self.gff_parser.all_seq))
            fout.write ("Valid sequence(s) in gff file\t{}\n".format(self.gff_parser.valid_seq))
            fout.write ("Features(s) in gff file\t{}\n".format(self.gff_parser.all_features))
            fout.write ("Valid features(s) in gff file\t{}\n".format(self.gff_parser.valid_features))
            if self.fusion:
                fout.write ("Remaining features after fusion\t{}\n".format(self.gff_parser.fused_features))

    ##### PRIVATE METHODS #####

    def extract_seq(self, seq_id, start, end, strand):

        if strand == "+":
            return str(self.seq_dict[seq_id][start+1:end+1].seq)
        else:
            return str(self.seq_dict[seq_id][start+1:end+1].reverse_complement().seq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    start_time = time()
    gf_extractor = GffFastaExtractor.class_init()
    gf_extractor()
    print ("\nDone in {}s".format(round(time()-start_time, 3)))
