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

# Third party package
from Bio import SeqIO

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GffFastaExtractor (object):

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "GffFastaExtractor 0.1"
    USAGE = "Usage: %prog -f FASTA_file -g GFF_file [...]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-f', '--fasta', dest="fasta",
            help= "Path to the fasta file contaning the complete genome")
        optparser.add_option('-g', '--gff', dest="gff",
            help= "Path to the gff file containing annotations of the genome file")
        optparser.add_option('-o', '--offset', dest="offset", default = 0,
            help= "Bases to extract before and after the feature (facultative, default 0)")
        optparser.add_option('--fusion', dest="fusion", action='store_true', default = False,
            help= "Fuse overlapping features in a meta-feature (facultative, default False)")
        optparser.add_option('--features', dest="features", default = "",
            help= "Restrict extraction to a list of features. The list must be SPACE separated and quoted(facultative, default all)")
        optparser.add_option('--chromosomes', dest="chromosomes", default = "",
            help= "Restrict extraction to a list of chromosomes. The list must be SPACE separated and quoted(facultative, default all)")

        ### Parse arguments
        opt, args = optparser.parse_args()

        # Verify the presence of options
        try:
            assert opt.fasta, "Missing fasta (-f) option"
            assert opt.gff, "Missing gff (-g) option"

        except AssertionError as E:
            print (E)
            optparser.print_help()
            sys_exit()

        ### Init a RefMasker object
        return GffFastaExtractor (opt.fasta, opt.gff, opt.offset, opt.fusion, opt.features.split(), opt.chromosomes.split())

        # optparser.print_help()

    ##### FONDAMENTAL METHODS #####

    def __init__ (self, fasta, gff, offset=0, fusion=False, features=[], chromosomes=[]):
        """Init fonction parsing fasta and gff files"""

        print("Initialize GffFastaExtractor")

        self.fasta = fasta
        self.gff = gff
        self.offset = int(offset)
        self.fusion = fusion
        self.features = features
        self.chromosomes = chromosomes
        self.separator = "|" # substitution separator

        self.out_name = "{}_{}_Offset-{}_Features-{}_Chr-{}_{}.fa.gz".format (
            self.fasta.rpartition('/')[2].partition('.')[0],
            self.gff.rpartition('/')[2].partition('.')[0],
            self.offset,
            "-".join(self.features) if self.features else "all",
            "-".join(self.chromosomes) if self.chromosomes else "all",
            "fused" if self.fusion else "not_fused")

        # Parse the fasta sequence and store in a dict
        print("  Parsing fasta file")
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
        print("  Parsing gff file")
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

        start_time = time()

        # Parse the gff parser and the sequence dictionary to extract the sequence of the features
        print("Extract features and write fasta output")

        n_feature = 0

        with gopen (self.out_name, "w") as fasta_out:
            for seq_id, gff_sequence in self.gff_parser.gff_dict.items():

                assert seq_id in self.seq_dict, "fasta and gff are incompatible: {} not found in fasta".format(seq_id)

                print ("  Extracting features from sequence {}".format(seq_id))
                for feature in gff_sequence.features:

                    fasta_out.write(">{}\n{}\n".format(
                            str(feature).replace("\t", self.separator).replace(" ", "_"),
                            self.extract_seq(seq_id, feature.start+1, feature.end, feature.strand)))
                    n_feature += 1

        # Finalize
        print ("\nExtract and wrote {} features in {}s".format(n_feature, round(time()-start_time, 3)))
        return(0)

    ##### PRIVATE METHODS #####

    def extract_seq(self, seq_id, start, end, strand):

        if strand == "+":
            return str(self.seq_dict[seq_id][start:end].seq)
        else:
            return str(self.seq_dict[seq_id][start:end].reverse_complement().seq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GffFeature(object):
    """Object representing a gff feature"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __init__(self, seq_id, source, type, start, end, score, strand, phase, attributes):
        self.seq_id = seq_id
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __str__ (self):
        """Return a gff formated line"""
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.seq_id,
            self.source,
            self.type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.phase,
            ";".join(["{}={}".format(key, ":".join(attribute)) for key, attribute in self.attributes.items()]))

    def __add__(self, other):
        """Support for concatenation of GffFeature objects = + operator"""
        assert self.seq_id == other.seq_id, "Sequence ID are not equal"
        assert self.type == other.type, "Feature types are not equal"
        assert self.strand == other.strand, "Feature strand are not equal"
        assert self.attributes["gene_id"] == other.attributes["gene_id"], "Feature Gene ID are not equal"

        # Fuse attribute dictionnary
        attributes_dict = OrderedDict()
        for key in self.attributes.keys():
            attributes_dict[key] = list(set(self.attributes[key]+other.attributes[key]))

        return GffFeature(
            seq_id = self.seq_id,
            source = self.source if self.source == other.source else "{}:{}".format (self.source, other.source),
            type = self.type,
            start = self.start if self.start <= other.start else other.start,
            end = self.end if self.end >= other.end else other.end,
            score = '.',
            strand = self.strand,
            phase = '.',
            attributes = attributes_dict)

    def is_overlapping(self, other):
        """Verify if 2 features object overlap or are adjacent"""
        assert self.attributes["gene_id"] == other.attributes["gene_id"], "Feature Gene ID are not equal"
        return self.seq_id == other.seq_id and self.start <= other.end+1 and other.start <= self.end+1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GffSequence(object):
    """Object representing a sequence in a gff file"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __init__(self, seq_id, start, end):
        self.seq_id = seq_id
        self.start = start
        self.end = end
        self.features = []
        self.features_count = 0

    def __str__ (self):
        """Return a gff formated line"""
        return "##sequence-region {} {} {}\n{}".format(
            self.seq_id,
            self.start,
            self.end,
            "\n".join([str(feature) for feature in self.features]))

    def add_feature(self, source, type, start, end, score, strand, phase, attributes):
        self.features.append(GffFeature(self.seq_id, source, type, start, end, score, strand, phase, attributes))
        self.features_count += 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GffParser(object):
    """ A minimalistic GFF3 format parser supporting gzip decompression """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    ##### FONDAMENTAL METHODS #####

    def __init__(self, gff_file, fusion=False, offset=0, features=[], chromosomes=[]):
        """
        gff_file    Path to a standard gff3 file
        features    restrict to the type of features indicated in the list (Default = all types)
        """
        openFunc = gopen if gff_file.endswith(".gz") else open
        self.gff_dict = OrderedDict()
        self.gff_header = ""
        self.all_features = self.valid_features = self.fused_features = self.all_seq = self.valid_seq = 0

        with openFunc(gff_file) as gff:
            for line in gff:
                line = line.lstrip()
                if not line:
                    continue

                # Parse sequence delimiter in gff file (ex: ##sequence-region chr1 1 248956422)
                if line.startswith("##sequence-region"):

                    # Part and verify that the number of fields is standard
                    parts = line.strip().split()
                    assert len(parts) == 4, "File format is not standard-compatible: Invalid sequence-region pragma"

                    # Create a GffSequence
                    seq_id = unquote(parts[1])
                    self.all_seq += 1
                    if not chromosomes or seq_id in chromosomes:
                        self.valid_seq += 1
                        self.gff_dict[seq_id]= GffSequence (
                            seq_id = seq_id,
                            start = int(parts[2]),
                            end = int(parts[3]))

                # Parse the general gff header
                elif line.startswith("#"):
                    self.gff_header+=line

                # Parse gff features
                else:

                    # Count all features
                    self.all_features += 1
                    if self.all_features % 100000 == 0:
                        print ("  {} features parsed".format(self.all_features))

                    # Part and verify that the number of fields is standard
                    parts = line.strip().split("\t")
                    assert len(parts) == 9, "File format is not standard-compatible : Invalid feature description line"

                    # Verify first if the feature type and chromosome is allowed
                    type = '.' if parts[2] == "." else unquote(parts[2])
                    if not features or type in features:
                        seq_id = '.' if parts[0] == "." else unquote(parts[0])
                        if not chromosomes or seq_id in chromosomes:

                            # Count valid features
                            self.valid_features += 1

                            # Verify that the corresponding sequence-region pragma is already in the dict
                            assert seq_id != ".", "File format is not standard-compatible : Missing seq_id of a feature"
                            assert seq_id in self.gff_dict, "File format is not standard-compatible : Missing or misplaced sequence-region pragma"

                            # Extract start and end coordinates, apply offset and verify that the values are valid
                            start = int(parts[3]) - offset
                            if start < 0:
                                start = 0

                            end = int(parts[4]) + offset
                            if end > self.gff_dict[seq_id].end:
                                end = self.gff_dict[seq_id].end

                            # Finally add the feature to the feature list corresponding to the seqid entry in gff_dict
                            self.gff_dict[seq_id].add_feature(
                                source = '.' if parts[1] == "." else unquote(parts[1]),
                                type = type,
                                start = start,
                                end = end,
                                score = '.' if parts[5] == '.' else float(parts[4]),
                                strand = '.' if parts[6] == '.' else unquote(parts[6]),
                                phase = '.' if parts[7] == '.' else int(parts[7]),
                                attributes = self._parse_attributes(parts[8]))

        # If the fusion of the features is required = additional processing needed
        if fusion:
            print ("Fuse overlapping features")

            for seq_id in self.gff_dict.keys():

                # sort the list by start coordinates
                self.gff_dict[seq_id].features.sort(key=lambda feature: feature.start)

                # init empty features and a list for the new features
                positive_feature = None
                negative_feature = None
                new_feature_list = []

                for feature in self.gff_dict[seq_id].features:

                    # Process features on the positive strand
                    if feature.strand == "+":
                        if not positive_feature:
                            positive_feature = feature
                        else:
                            if feature.is_overlapping (positive_feature):
                                positive_feature = positive_feature+feature
                            else:
                                new_feature_list.append (positive_feature)
                                positive_feature = feature

                    # Process features on the negative strand
                    else :
                        if not negative_feature:
                            negative_feature = feature
                        else:
                            if feature.is_overlapping (negative_feature):
                                negative_feature = negative_feature+feature
                            else:
                                new_feature_list.append (negative_feature)
                                negative_feature = feature

                # Flush the last item if needed
                if positive_feature:
                    new_feature_list.append (positive_feature)
                if negative_feature:
                    new_feature_list.append (negative_feature)

                # Replace the previous feature list with the new one
                self.gff_dict[seq_id].features = new_feature_list
                self.fused_features += len (new_feature_list)


    def __str__ (self):
        """Return a gff formated line"""

        return "{}{}".format(
            self.gff_header,
            "\n".join([str(sequence) for sequence in self.gff_dict.values()]))

    ##### PRIVATE METHODS #####

    def _parse_attributes(self, attribute_string):
        """Parse the GFF3 attribute column and return a dict"""

        # Empty dict to collect features attributes
        attributes_dict = OrderedDict()

        # Parse and return a dict of all attributes if attribute_string is defined
        if attribute_string != ".":
            for attribute in attribute_string.split(";"):
                key, value = attribute.split("=")
                attributes_dict[unquote(key)] = unquote(value).split(" ") #split in a list of value
        return attributes_dict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    gf_extractor = GffFastaExtractor.class_init()
    gf_extractor()
