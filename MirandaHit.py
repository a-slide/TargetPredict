# -*- coding: utf-8 -*-

"""
@package    Miranda
@brief      Class to represent the informations from a miranda hit
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
from collections import OrderedDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class MirandaHit(object):
    """
    @class  MirandaHit
    @brief  Object oriented class containing informations of one miranda hit
    """
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
        """
        Create a MirandaHit object which is automatically added to the class tracking instance list
        """

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

    def __str__(self):
        msg = "HIT {}".format(self.id)
        msg += "\tQuery\t{}:{}-{}\n".format(self.q_id, self.q_start, self.q_end)
        msg += "\tSubject\t{}:{}-{}\n".format(self.s_id, self.s_start, self.s_end)
        msg += "\tLenght : {}\tIdentity : {}%\Homology : {}%\tScore : {}\tEnergy : {}\n".format(
            self.length, self.identity, self.homology, self.score, self.energy )
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def get_report (self, full=False):
        """
        Generate a report under the form of an Ordered dictionary
        @param full If true a dict containing all self parameters will be returned
        """
        report = OrderedDict ()
        report["Query"] = "{}:{}-{}".format(self.q_id, self.q_start, self.q_end)
        report["Subject"] = "{}:{}-{}".format(self.s_id, self.s_start, self.s_end)

        if full:
            report["Identity"] = self.identity
            report["Homology"] = self.homology
            report["Score"] = self.score
            report["Energy"] = self.energy
            report["Hit length"] = self.length

        return report
