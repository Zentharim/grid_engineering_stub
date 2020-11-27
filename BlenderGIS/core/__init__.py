# KEEP io_import_shp

import logging
logging.basicConfig(level=logging.getLevelName('INFO'))

from .checkdeps import HAS_GDAL, HAS_PYPROJ, HAS_IMGIO, HAS_PIL
from .settings import getSettings, setSettings
from .errors import OverlapError

from .utils import XY, BBOX

from .proj import SRS, Reproj, reprojPt, reprojPts, reprojBbox, reprojImg

from .lib import shapefile
