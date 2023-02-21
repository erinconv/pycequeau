from __future__ import annotations

import os
import numpy as np
import xarray as xr

_default_data_path = "pycequeau/test/data/"
_default_cache_dir = "titorial_data"

_internal_files =dict(
    precipitation="pr_NAM-44_CCCma-CanESM2_rcp85_r1i1p1_UQAM-CRCM5_v1_day.nc",
)


def open_example(name:str):
    if name in _internal_files:        
        file_path = os.path.join(_default_data_path,
                                 _internal_files[name])
    else:
        raise ValueError("The example file does not exist")
    
    # Check the file extension
    if file_path.endswith(".nc"):
        ds = xr.open_dataset(file_path)
        return ds["pr"]
