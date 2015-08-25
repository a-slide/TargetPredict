# -*- coding: utf-8 -*-

"""
@package    BlastHit
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library packages
from collections import OrderedDict

# Local packages
from GffLine import GffLine

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BlastHit(object):
    """ Object oriented class containing informations of one blast hit """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#
    ID_COUNT = 0

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def next_id (self):
        cur_id = self.ID_COUNT
        self.ID_COUNT +=1
        return cur_id

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, q_id, s_id, identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, score, q_seq):
        """ Create a BlastHit object """

        # Store parameters in self variables
        self.id = self.next_id()
        self.q_id = q_id
        self.s_id = s_id
        self.identity = float(identity)
        self.length = int(length)
        self.mis = int(mis)
        self.gap = int(gap)
        self.evalue = float(evalue)
        self.score = float(score)
        self.q_seq = q_seq
        self.q_start = int(q_start)-1
        self.q_end = int(q_end)

        # Parse the gff_line in the contained in the subject id
        self.gff = GffLine(self.s_id)

        # Correct coordinates of hit for python 0 based coordinates depending of the orientation

        if int(s_start) < int(s_end):
            self.strand = "+"
            self.s_start = int(s_start)-1
            self.s_end = int(s_end)
        else:
            self.strand = "-"
            self.s_start = int(s_end)-1
            self.s_end = int(s_start)

    def __str__(self):
        msg = "HIT {}".format(self.id)
        msg += "\tQuery\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.q_end)
        msg += "\tSubject\t{}:{}-{}({})\n".format(self.s_id, self.s_start, self.s_end, self.strand)
        msg += "\tLenght : {}\tIdentity : {}%\tEvalue : {}\tBit score : {}\n".format(self.length, self.identity, self.evalue, self.score)
        msg += "\tAligned query seq : {}\n".format(self.q_seq)
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def report (self):
        """Return and ordered dict containing all the required informations"""
        report = OrderedDict()
        for key, val in self.gff.report().items():
            report ["subject_"+key] = val
        report ["subject_hit_start"] = self.s_start+1
        report ["subject_hit_end"] = self.s_end
        report ["query"] = self.q_id
        report ["query_hit_start"] = self.q_start+1
        report ["query_hit_end"] = self.q_end
        report ["hit_score"] = self.score
        report ["hit_evalue"] = self.evalue
        report ["hit_length"] = self.length
        report ["hit_identity"] = self.identity
        report ["hit_gap"] = self.gap
        report ["hit_mismatch"] = self.mis

        return report
