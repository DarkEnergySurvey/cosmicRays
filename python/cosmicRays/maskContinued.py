
#
# LSST Data Management System
# Copyright 2008-2017 LSST/AURA.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ["Mask", "MaskPixel"]

import numpy as np

from cosmicRays.utils import TemplateMeta
from cosmicRays.fitsIoWithOptions import imageReadFitsWithOptions, imageWriteFitsWithOptions
from libcosmicRays.afw.image import MaskX
from cosmicRays.slicing import supportSlicing
from cosmicRays.disableArithmetic import disableMaskArithmetic

MaskPixel = np.int32


class Mask(metaclass=TemplateMeta):
    TEMPLATE_PARAMS = ("dtype",)
    TEMPLATE_DEFAULTS = (MaskPixel,)

    def __reduce__(self):
        from lsst.afw.fits import reduceToFits
        return reduceToFits(self)

    def __str__(self):
        return "{}, bbox={}, maskPlaneDict={}".format(self.array, self.getBBox(), self.getMaskPlaneDict())

    def __repr__(self):
        return "{}.{}={}".format(self.__module__, self.__class__.__name__, str(self))

    readFitsWithOptions = classmethod(imageReadFitsWithOptions)

    writeFitsWithOptions = imageWriteFitsWithOptions


Mask.register(MaskPixel, MaskX)
Mask.alias("X", MaskX)

for cls in (MaskX, ):
    supportSlicing(cls)
    disableMaskArithmetic(cls)
