# -*- coding: utf-8 -*-

"""
@package RefMasker
@brief  Collection of helper functions for file and path manipulation
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
from os import access, R_OK, path
from gzip import open as gopen
from shutil import copy

#~~~~~~~ PREDICATES ~~~~~~~#

def is_readable_file (fp):
    """ Verify the readability of a file or list of file """
    return access(fp, R_OK)

def is_gziped (fp):
    """ Return True if the file is Gziped else False """
    return fp[-2:].lower() == "gz"

#~~~~~~~ PATH MANIPULATION ~~~~~~~#

def file_basepath (fp):
    """Return the path of a file without the last extension """
    return fp.rpartition('.')[0]

def file_basename (fp):
    """Return the basename of a file without folder location and extension """
    return fp.rpartition('/')[2].partition('.')[0]

def file_extension (fp):
    """ Return The extension of a file in lowercase """
    return fp.rpartition(".")[2].lower()

def file_name (fp):
    """ Return The complete name of a file with the extension but without folder location """
    return fp.rpartition("/")[2]

def dir_name (fp):
    """ Return the complete path where is located the file without the file name """
    return fp.rpartition("/")[0].rpartition("/")[2]

def rm_blank (name, replace=""):
    """ Replace blank spaces in a name by a given character (default = remove)
    Blanks at extremities are always removed and nor replaced """
    return replace.join(name.split())

#~~~~~~~ FILE MANIPULATION ~~~~~~~#

def gunzip (src, dst):
    """
    @param source Path of the input compressed file
    @param destination Path to a directory or a file where source will be extracted. If destination
    is a directory, the file will be extracted into and automatically renamed.
    @return The path of the destination file
    """
    # Generate a automatic name without .gz extension in the directory
    if path.isdir(dst):
        dst = path.join(dst, file_name(src.rstrip(".gz")))

    # Try to initialize handle for the compressed file
    with gopen(src, 'rb') as in_handle:
        with open(dst, "wb") as out_handle:
        # Write input file in output file
            out_handle.write (in_handle.read())

    return dst

def cp (src, dst):
    """
    @param source Path of the input file
    @param destination Path to a directory or a file where source will be copied. If destination
    is a directory, the file will be copied into and automatically renamed.
    @return The path of the destination file
    """
    # Generate a automatic name in the directory from the source file name
    if path.isdir(dst):
        dst = path.join(dst, file_name(src))

    copy(src, dst)

    return dst
