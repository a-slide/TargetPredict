# -*- coding: utf-8 -*-

"""
@package    GffLine
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library packages
from collections import OrderedDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GffLine(object):
    """
    @class  GffLine
    @brief  Parse and store informations of one GffLine
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, gff_line):
        """
        Create a MirandaHit object which is automatically added to the class tracking instance list
        """

        self.separator = "#" # substitution separator

        # Store parameters in self variables
        parts = gff_line.strip().split(self.separator)
        assert len(parts) == 9, "File format is not standard-compatible:\n{}\n{}".format(gff_line, parts)

        self.seq_id = parts[0]
        self.source = parts[1]
        self.type = parts[2]
        self.start = parts[3]
        self.end = parts[4]
        self.score = parts[5]
        self.strand = parts[6]
        self.phase = parts[7]

        # Add attibutes names in the object self dict
        for attribute in parts[8].split(";"):
            key, value = attribute.split("=")
            self.__dict__[key] = value

    def __str__(self):
        return "".join (["{}\t{}\n".format (key, value) for key, value in self.__dict__.items ()])

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def report (self):
        report = OrderedDict()
        report ["Exon ID"] = self.exon_id
        report ["Exon Number"] = self.exon_number
        report ["Transcript ID"] = self.transcript_id
        report ["Transcript Name"] = self.transcript_name
        report ["Gene ID"] = self.gene_id
        report ["Gene Name"] = self.gene_name
        report ["Chromosome"] = self.seq_id
        report ["Start"] = self.start
        report ["End"] = self.end
        report ["Strand"] = self.strand

        return report

# chrX#HAVANA#exon#31507281#31507453#.#-#.#ID=exon:ENST00000358062.6:9;Parent=ENST00000358062.6;gene_id=ENSG00000198947.14;transcript_id=ENST00000358062.6;gene_type=protein_coding;gene_status=KNOWN;gene_name=DMD;transcript_type=protein_coding;transcript_status=KNOWN;transcript_name=DMD-017;exon_number=9;exon_id=ENSE00003500236.1;level=2;protein_id=ENSP00000350765.2;transcript_support_level=5;havana_gene=OTTHUMG00000021336.5;havana_transcript=OTTHUMT00000355787.1;tag=mRNA_start_NF,cds_start_NF
