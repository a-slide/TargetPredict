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
                            start = int(parts[3])-1
                            end = int(parts[4])-1
                            attributes = self._parse_attributes(parts[8])

                            if offset:
                                start -= offset
                                end += offset
                                if start < 0:
                                    start = 0
                                if end > self.gff_dict[seq_id].end-1:
                                    end = self.gff_dict[seq_id].end-1
                                attributes["offset"] = [str(offset)]

                            # Finally add the feature to the feature list corresponding to the seqid entry in gff_dict
                            self.gff_dict[seq_id].add_feature(
                                source = '.' if parts[1] == "." else unquote(parts[1]),
                                type = type,
                                start = start,
                                end = end,
                                score = '.' if parts[5] == '.' else float(parts[4]),
                                strand = '.' if parts[6] == '.' else unquote(parts[6]),
                                phase = '.' if parts[7] == '.' else int(parts[7]),
                                attributes = attributes)

        # If the fusion of the features is required = additional processing needed
        if fusion:
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
class GffFeature(object):
    """Object representing a gff feature"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __init__(self, seq_id, source, type, start, end, score, strand, phase, attributes, fused_features=1):
        self.seq_id = seq_id
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self.fused_features = fused_features

    @property
    def attribute_string (self):
        return ";".join(["{}={}".format(key, ",".join(attribute)) for key, attribute in self.attributes.items()])

    def __str__ (self):
        """Return a gff formated line"""

        if self.fused_features >= 2 :
            return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{};fused_features={}".format (self.seq_id,
                self.source, self.type, self.start+1, self.end+1, self.score, self.strand, self.phase,
                self.attribute_string, self.fused_features)
        else:
            return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format (self.seq_id, self.source, self.type,
            self.start+1, self.end+1, self.score, self.strand, self.phase, self.attribute_string)

    def __add__(self, other):
        """Support for concatenation of GffFeature objects = + operator"""
        assert self.seq_id == other.seq_id, "Sequence ID are not equal: {} - {}".format(
            self.seq_id, other.seq_id)
        assert self.type == other.type, "Feature types are not equal: {} - {}".format(
            self.type, other.type)
        assert self.strand == other.strand, "Feature strand are not equal: {} - {}".format(
            self.strand, other.strand)

        # Fuse attribute dictionnaries
        attributes_dict = OrderedDict()
        # start by the keys of self
        for key in self.attributes.keys():
            if key not in other.attributes:
                attributes_dict[key] = self.attributes[key]
            else:
                attributes_dict[key] = list(set(self.attributes[key]+other.attributes[key]))
        # Add the orphan keys of other
        for key in other.attributes.keys():
            if key not in self.attributes:
                attributes_dict[key] = other.attributes[key]

        # finally return the fused GffFeature object
        return GffFeature(
            seq_id = self.seq_id, source = self.source, type = self.type,
            start = self.start if self.start <= other.start else other.start,
            end = self.end if self.end >= other.end else other.end,
            score = '.', strand = self.strand, phase = '.', attributes = attributes_dict,
            fused_features = self.fused_features + other.fused_features)

    def is_overlapping(self, other):
        """Verify if 2 features object overlap or are adjacent"""
        return self.seq_id == other.seq_id and self.start <= other.end and other.start <= self.end
