from __future__ import annotations

import numpy as np
from osgeo import gdal,ogr


class Basin:
    def __init__(self,
                 name: str,
                 DEM: gdal.Dataset,
                 FAC: gdal.Dataset,
                 CEgrid: gdal.Dataset,
                 fishnet: ogr.DataSource,
                 watershed: ogr.DataSource,
                 grid_size: int) -> None:
        self.watershed = watershed
        self.fishnet = fishnet
        self.DEM = DEM
        self.CEgrid = CEgrid
        self.grid_size = grid_size