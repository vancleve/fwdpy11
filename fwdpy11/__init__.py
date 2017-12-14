#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
#

import sys

if sys.version_info[0] < 3:
    raise ValueError("Python3 required!")

from fwdpy11._version import __version__ # NOQA

from .fwdpp_types import *
from ._opaque_gametes import *
from ._opaque_mutations import *
from ._opaque_diploids import *
from ._opaque_generalmutvecs import *
from .fwdpy11_types import *
from ._regions import *
from ._dev import *
from ._gslrng import GSLrng

class SingleLocusDiploid(fwdpy11_types.SingleLocusDiploid):
    """
    Diploid data type for a single (usually contiguous) genomic region
    """
    pass

class SlocusPop(fwdpy11_types.SlocusPop):
    """
    Population object representing a single 
    deme and a single genomic region.
    """
    pass


class MlocusPop(fwdpy11_types.MlocusPop):
    """
    Representation of a multi-locus, single 
    deme system.
    """
    pass


class SlocusPopGeneralMutVec(fwdpy11_types.SlocusPopGeneralMutVec):
    """
    Single-deme object using 
    :class:`fwpy11.GeneralMutVec`
    as the mutation type.
    """
    pass
