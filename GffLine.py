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

        self.separator = "|" # substitution separator

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
        report ["exon_id"] = self.exon_id
        report ["exon_number"] = self.exon_number
        report ["transcript_id"] = self.transcript_id
        report ["transcript_name"] = self.transcript_name
        report ["gene_id"] = self.gene_id
        report ["gene_name"] = self.gene_name
        report ["chromosome"] = self.seq_id
        report ["start"] = self.start
        report ["end"] = self.end
        report ["strand"] = self.strand

        return report
