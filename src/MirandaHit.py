# -*- coding: utf-8 -*-

"""
@package    MirandaHit
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
class MirandaHit(object):
    """ Object oriented class containing informations of one miranda hit """
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

    def __init__(self, q_id, s_id, score, energy, q_start, q_end, s_start, s_end, length, identity, homology):
        """ Create a MirandaHit object """

        # Store parameters in self variables
        self.id = self.next_id()
        self.q_id = q_id
        self.s_id = s_id
        self.score = float(score)
        self.energy = float(energy)
        self.q_start = int(q_start)
        self.q_end = int(q_end)
        self.s_start = int(s_start)
        self.s_end = int(s_end)
        self.length = int(length)
        self.identity = float(identity[:-1])
        self.homology = float(homology[:-1])

        # Parse the gff_line in the contained in the subject id
        self.gff = GffLine(self.s_id)

    def __str__(self):
        msg = "HIT {}".format(self.id)
        msg += "\tQuery\t{}:{}-{}\n".format(self.q_id, self.q_start, self.q_end)
        msg += "\tSubject\t{}:{}-{}\n".format(self.s_id, self.s_start, self.s_end)
        msg += "\tLenght : {}\tIdentity : {}%\Homology : {}%\tScore : {}\tEnergy : {}\n".format(
            self.length, self.identity, self.homology, self.score, self.energy )
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
        report ["hit_energy"] = self.energy
        report ["hit_length"] = self.length
        report ["hit_identity"] = self.identity
        report ["hit_homology"] = self.homology

        return report
