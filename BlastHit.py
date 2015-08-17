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

        # Parse the gff_line in the contained in the subject id
        self.gff = GffLine(self.s_id)

        # Correct coordinates of hit for python 0 based coordinates depending of the orientation
        if int(q_start) < int(q_end):
            self.q_orient = "+"
            self.q_start = int(q_start)-1
            self.q_end = int(q_end)
        else:
            self.q_orient = "-"
            self.q_start = int(q_start)
            self.q_end = int(q_end)-1

        if int(s_start) < int(s_end):
            self.s_orient = "+"
            self.s_start = int(s_start)-1
            self.s_end = int(s_end)
        else:
            self.s_orient = "-"
            self.s_start = int(s_start)
            self.s_end = int(s_end)-1

    def __str__(self):
        msg = "HIT {}".format(self.id)
        msg += "\tQuery\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.q_end, self.q_orient)
        msg += "\tSubject\t{}:{}-{}({})\n".format(self.s_id, self.s_start, self.s_end, self.s_orient)
        msg += "\tLenght : {}\tIdentity : {}%\tEvalue : {}\tBit score : {}\n".format(self.length, self.identity, self.evalue, self.score)
        msg += "\tAligned query seq : {}\n".format(self.q_seq)
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def report (self):
        report = OrderedDict()
        report ["Query"] = self.q_id
        report ["Query_start"] = self.q_start
        report ["Query_end"] = self.q_end
        report ["Subject"] = self.gff.report()
        report ["Subject_start"] = self.s_start
        report ["Subject_end"] = self.s_end
        report ["Score"] = self.score
        report ["Evalue"] = self.evalue
        report ["Length"] = self.length
        report ["Identity"] = self.identity

        return report
