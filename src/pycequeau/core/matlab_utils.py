"""Utilities for converting between Python data and MATLAB (.mat) formats."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.io.matlab import mat_struct



def to_mat_compatible(value):
    """Recursively convert Python containers to scipy.savemat-compatible values."""
    if isinstance(value, dict):
        return {k: to_mat_compatible(v) for k, v in value.items()}
    if isinstance(value, (set, tuple)):
        value = list(value)
    if isinstance(value, list):
        if not value:
            return np.array([])
        if all(isinstance(v, (int, float, np.integer, np.floating, bool)) for v in value):
            return np.asarray(value)
        return np.array([to_mat_compatible(v) for v in value], dtype=object)
    return value


def dataframe_to_struct_array(df: pd.DataFrame) -> np.ndarray:
    """Convert DataFrame to MATLAB struct-array layout (1xN struct)."""
    field_names = [str(col) for col in df.columns]
    struct_dtype = np.dtype([(name, object) for name in field_names])
    struct_arr = np.empty((1, len(df)), dtype=struct_dtype)

    for row_idx, (_, row) in enumerate(df.iterrows()):
        for field in field_names:
            struct_arr[field][0, row_idx] = mat_field_value(row[field])

    return struct_arr


def mat_field_value(value):
    """Normalize a single table value to MATLAB-friendly field content."""
    if isinstance(value, (set, tuple)):
        value = list(value)
    if isinstance(value, (np.ndarray, pd.Series)):
        arr = np.asarray(value)
        if arr.dtype == object:
            return np.array([mat_field_value(v) for v in arr], dtype=object)
        return arr
    if isinstance(value, list):
        if all(isinstance(v, (int, float, np.integer, np.floating, bool)) for v in value):
            return np.asarray(value)
        return np.array([mat_field_value(v) for v in value], dtype=object)
    if isinstance(value, np.generic):
        return value.item()
    return value


def mat_to_py(value):
    """Recursively convert scipy.loadmat objects into plain Python containers."""
    # scipy.io MATLAB struct: use vars() (public __dict__) instead of _fieldnames
    if isinstance(value, mat_struct):
        return {
            k: mat_to_py(v)
            for k, v in vars(value).items()
            if not k.startswith("_")
        }
    if isinstance(value, np.ndarray):
        if value.dtype == object:
            return [mat_to_py(v) for v in value.tolist()]
        return value.tolist()
    return value
