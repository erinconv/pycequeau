from __future__ import annotations

import os
from dataclasses import dataclass

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from osgeo import gdal, ogr
from shapely.geometry import Point

from ..core import manage_files, projections, utils as u
from ..physiographic.base import Basin
from .base import Meteo

__methods__ = [
    "linear",
    "nearest",
]


@dataclass(frozen=True)
class NetCDFGridConfig:
    """Configuration for the NetCDF-to-CEQUEAU workflow."""

    time_name: str = "time"
    lat_name: str = "lat"
    lon_name: str = "lon"
    ce_index_offset: int = 10


class StationNetCDF(Meteo):
    """Meteorological workflow for gridded NetCDF datasets."""

    def __init__(
        self,
        basin_struct: Basin,
        ds: xr.Dataset,
        config: NetCDFGridConfig | None = None,
    ) -> None:
        super().__init__(basin_struct)
        self.config = config or NetCDFGridConfig()
        self.ds = _standardize_dataset(ds, self.config)
        self.table: pd.DataFrame | None = None
        self.lon_utm: np.ndarray | None = None
        self.lat_utm: np.ndarray | None = None
        self.interpolated: xr.Dataset | None = None

        ce_area = float(self.basin_struct.bassinVersant["superficieCE"]) * 1e6
        self.basin_struct.set_dimensions(np.sqrt(ce_area), np.sqrt(ce_area))

    @classmethod
    def from_dataset(
        cls,
        basin_struct: Basin,
        ds: xr.Dataset,
        config: NetCDFGridConfig | None = None,
    ) -> "StationNetCDF":
        return cls(basin_struct, ds, config=config)

    @classmethod
    def from_netcdf_folder(
        cls,
        basin_struct: Basin,
        vars_path: str,
        source: str = "ERA",
        config: NetCDFGridConfig | None = None,
    ) -> "StationNetCDF":
        source_name = source.upper()
        vars_dict = manage_files.dict_netCDF(vars_path)
        if source_name == "ERA":
            ds = manage_files.get_ERA_Dataset(vars_dict)
        elif source_name == "CORDEX":
            ds = manage_files.get_CORDEX_Dataset(vars_dict)
        else:
            raise ValueError(
                f"Unsupported NetCDF source '{source}'. Expected 'ERA' or 'CORDEX'."
            )
        return cls(basin_struct, ds, config=config)

    @classmethod
    def charge_ERA_Meteo(cls, bassinVersant: Basin, vars_path: str) -> "StationNetCDF":
        """Backward-compatible constructor for ERA datasets."""
        return cls.from_netcdf_folder(bassinVersant, vars_path, source="ERA")

    @classmethod
    def charge_CORDEX_Meteo(
        cls, bassinVersant: Basin, vars_path: str
    ) -> "StationNetCDF":
        """Backward-compatible constructor for CORDEX datasets."""
        return cls.from_netcdf_folder(bassinVersant, vars_path, source="CORDEX")

    def stations_table(
        self,
        name_meteo_grid_file: str = "meteo_grid_points.shp",
        export: bool = True,
    ) -> pd.DataFrame:
        """Build the station table used for raster interpolation."""
        dem = gdal.Open(self.basin_struct.DEM, gdal.GA_ReadOnly)
        ce_grid = gdal.Open(self.basin_struct.get_CEgrid, gdal.GA_ReadOnly)
        watershed = ogr.Open(self.basin_struct.watershed_shapefile, gdal.GA_ReadOnly)
        epsg_dem = projections.get_proj_code(dem)

        xy_pair = _get_netcdf_grid_points(self.ds, ce_grid, watershed, self.config)
        lon_utm, lat_utm = projections.latlon_to_utm(
            xy_pair[:, 1],
            xy_pair[:, 0],
            epsg_dem,
        )
        self.lon_utm = lon_utm
        self.lat_utm = lat_utm
        self.table = _create_station_table(
            ce_grid,
            dem,
            lon_utm,
            lat_utm,
            xy_pair,
            self.config,
        )

        if export:
            self._export_station_table(name_meteo_grid_file, epsg_dem)
        return self.table

    def interpolation(
        self,
        method: str = "nearest",
        name_meteo_grid_file: str = "meteo_grid_points.shp",
    ) -> xr.Dataset:
        """Interpolate the input dataset over the CE grid."""
        if method not in __methods__:
            raise ValueError(
                f"Unsupported interpolation method '{method}'. Expected one of {__methods__}."
            )
        if self.table is None:
            self.stations_table(name_meteo_grid_file=name_meteo_grid_file)

        self.interpolated = _interpolate_dataset_to_ce_grid(
            self.ds,
            self.basin_struct.get_CEgrid,
            self.table,
            method,
            self.config,
        )
        return self.interpolated

    def to_cequeau_grid(self, ds: xr.Dataset | None = None) -> xr.Dataset:
        """Export an interpolated dataset to the CEQUEAU meteorological layout."""
        source = ds if ds is not None else self.interpolated
        if source is None:
            raise ValueError("No interpolated dataset is available to convert.")
        return self.cequeau_grid(source, self.basin_struct)

    @classmethod
    def _cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        ce_df = _load_ce_table(basin_struct)
        return _build_cequeau_grid(ds, ce_df)

    def _export_station_table(self, file_name: str, epsg_code: int) -> None:
        if self.table is None:
            raise ValueError("The station table must be created before exporting it.")
        points = [
            Point(x, y)
            for y, x in zip(self.table["lat_utm"], self.table["lon_utm"])
        ]
        gdf = gpd.GeoDataFrame(self.table.copy(), geometry=points, crs=epsg_code)
        output_path = os.path.join(
            self.basin_struct.project_path,
            "geographic",
            file_name,
        )
        gdf.to_file(output_path)


def _standardize_dataset(ds: xr.Dataset, config: NetCDFGridConfig) -> xr.Dataset:
    rename_map: dict[str, str] = {}
    if "longitude" in ds.coords and config.lon_name not in ds.coords:
        rename_map["longitude"] = config.lon_name
    if "latitude" in ds.coords and config.lat_name not in ds.coords:
        rename_map["latitude"] = config.lat_name
    if rename_map:
        ds = ds.rename(rename_map)

    required_coords = {config.time_name, config.lat_name, config.lon_name}
    missing = required_coords.difference(ds.coords)
    if missing:
        missing_names = ", ".join(sorted(missing))
        raise ValueError(f"Missing NetCDF coordinates: {missing_names}")

    ds = ds.sortby(config.time_name)
    ds = ds.sortby(config.lat_name)
    ds = ds.sortby(config.lon_name)
    return ds


def _load_ce_table(basin_struct: Basin) -> pd.DataFrame:
    ce_path = os.path.join(
        basin_struct.project_path,
        "results",
        "carreauxEntiers.csv",
    )
    ce_df = pd.read_csv(ce_path, index_col=0)
    ce_df.index = ce_df["CEid"].astype(int).values
    ce_df["i"] = ce_df["i"].astype(int)
    ce_df["j"] = ce_df["j"].astype(int)
    return ce_df


def _build_cequeau_grid(ds: xr.Dataset, ce_df: pd.DataFrame) -> xr.Dataset:
    variable_names = [name for name in ds.data_vars if name != "CE"]
    if not variable_names:
        raise ValueError("The interpolated dataset does not contain meteorological variables.")

    time_values = pd.to_datetime(ds["time"].values)
    datenum = np.array(
        [366.0 + timestamp.toordinal() for timestamp in time_values],
        dtype=np.float32,
    )

    result = xr.Dataset(coords={
        "CEid": ce_df["CEid"].values.astype(np.int32),
        "pasTemp": datenum,
    })
    for variable_name in variable_names:
        result[variable_name] = (
            ("pasTemp", "CEid"),
            _extract_ce_timeseries(ds, ce_df, variable_name),
        )
        result[variable_name].attrs = ds[variable_name].attrs
    return result


def _extract_ce_timeseries(
    ds: xr.Dataset,
    ce_df: pd.DataFrame,
    variable_name: str,
) -> np.ndarray:
    ce_data = np.zeros((ds.sizes["time"], len(ce_df)), dtype=np.float32)
    for position, (_, ce_row) in enumerate(ce_df.iterrows()):
        i_index = int(ce_row["i"])
        j_index = int(ce_row["j"])
        ce_data[:, position] = ds[variable_name].loc[:, j_index, i_index].values.astype(
            np.float32
        )
    return ce_data


def _interpolate_dataset_to_ce_grid(
    ds: xr.Dataset,
    ce_grid_path: str,
    table: pd.DataFrame,
    method: str,
    config: NetCDFGridConfig,
) -> xr.Dataset:
    ce_grid = gdal.Open(ce_grid_path, gdal.GA_ReadOnly)
    ce_array = np.array(ce_grid.GetRasterBand(1).ReadAsArray())
    row_count, col_count = ce_array.shape
    i = np.arange(col_count, dtype=np.int16) + config.ce_index_offset
    j = np.arange(row_count, dtype=np.int16) + config.ce_index_offset

    working_ds = ds.copy()
    if float(working_ds[config.lon_name].max()) > 180:
        working_ds = working_ds.assign_coords(
            {config.lon_name: working_ds[config.lon_name].values - 360}
        )
        working_ds = working_ds.sortby(config.lon_name)

    lon_max_idx = u.find_nearest(working_ds[config.lon_name].values, table["lon"].max())
    lon_min_idx = u.find_nearest(working_ds[config.lon_name].values, table["lon"].min())
    lat_max_idx = u.find_nearest(working_ds[config.lat_name].values, table["lat"].max())
    lat_min_idx = u.find_nearest(working_ds[config.lat_name].values, table["lat"].min())

    working_ds = working_ds.isel(
        {
            config.lat_name: slice(min(lat_min_idx, lat_max_idx), max(lat_min_idx, lat_max_idx)),
            config.lon_name: slice(min(lon_min_idx, lon_max_idx), max(lon_min_idx, lon_max_idx)),
        }
    )

    working_ds = working_ds.assign_coords(
        {
            config.lon_name: np.linspace(i.min(), i.max(), working_ds.sizes[config.lon_name]),
            config.lat_name: np.linspace(j.min(), j.max(), working_ds.sizes[config.lat_name]),
        }
    )

    flipped_j = np.flip(j)
    interpolated = working_ds.interp(
        {
            config.time_name: working_ds[config.time_name],
            config.lat_name: flipped_j,
            config.lon_name: i,
        },
        method=method,
    )
    interpolated = interpolated.where(ce_array > 0)
    interpolated = interpolated.rename(
        {
            config.time_name: "time",
            config.lat_name: "j",
            config.lon_name: "i",
        }
    )
    interpolated = interpolated.transpose("time", "j", "i")
    interpolated["i"] = i.astype(np.int16)
    interpolated["j"] = flipped_j.astype(np.int16)
    interpolated = _append_ce_grid(interpolated, ce_grid)
    return interpolated.assign_attrs(
        interpolated=f"Interpolated using xarray.Dataset.interp with method='{method}'"
    )


def _get_netcdf_grid_points(
    ds: xr.Dataset,
    ce_grid: gdal.Dataset,
    watershed: ogr.DataSource,
    config: NetCDFGridConfig,
) -> np.ndarray:
    xtup, ytup, _ = u.GetExtent(ce_grid)
    epsg_dem = projections.get_proj_code(ce_grid)
    x, y = projections.utm_to_latlon(
        (np.amin(xtup), np.amax(xtup)),
        (np.amin(ytup), np.amax(ytup)),
        epsg_dem,
    )

    dy = abs(ds[config.lat_name][0].values - ds[config.lat_name][1].values)
    dx = abs(ds[config.lon_name][0].values - ds[config.lon_name][1].values)
    lon_mask = np.arange(x[0], x[1], dx)
    lat_mask = np.arange(y[0], y[1], dy)

    watershed_layer = watershed.GetLayer()
    xmin, xmax, ymin, ymax = watershed_layer.GetExtent()
    x, y = projections.utm_to_latlon((xmin, xmax), (ymin, ymax), epsg_dem)
    watershed_extent = (x, y)
    return u.falls_in_extent(watershed_extent, lon_mask, lat_mask)


def _create_station_table(
    ce_grid: gdal.Dataset,
    dem: gdal.Dataset,
    lon_utm: np.ndarray,
    lat_utm: np.ndarray,
    xy_pair: np.ndarray,
    config: NetCDFGridConfig,
) -> pd.DataFrame:
    ce_array = ce_grid.ReadAsArray()
    row_count, col_count = ce_array.shape
    i = np.arange(col_count, dtype=int)
    j = np.flip(np.arange(row_count, dtype=int))

    row, col = u.get_index_list(ce_grid, lon_utm, lat_utm)
    row = np.array(row, dtype=np.float32)
    col = np.array(col, dtype=np.float32)

    valid_row = row < row_count
    xy_pair = xy_pair[valid_row, :]
    col = col[valid_row]
    lat_utm = lat_utm[valid_row]
    lon_utm = lon_utm[valid_row]
    row = row[valid_row]

    valid_col = col < col_count
    xy_pair = xy_pair[valid_col, :]
    row = row[valid_col]
    lat_utm = lat_utm[valid_col]
    lon_utm = lon_utm[valid_col]
    col = col[valid_col]

    row = row.astype(np.int16)
    col = col.astype(np.int16)

    return pd.DataFrame(
        data={
            "id": [f"NC-grid-{num}" for num in range(len(i[col]))],
            "i": i[col] + config.ce_index_offset,
            "j": j[row] + config.ce_index_offset,
            "lat": xy_pair[:, 1],
            "lon": xy_pair[:, 0],
            "lat_utm": lat_utm,
            "lon_utm": lon_utm,
            "CEid": ce_array[row, col],
            "altitude": u.get_altitude_point(dem, lat_utm, lon_utm),
        }
    )


def _append_ce_grid(ds: xr.Dataset, ce_grid: gdal.Dataset) -> xr.Dataset:
    grid = ce_grid.ReadAsArray().astype(np.float16)
    grid[grid == 0] = np.nan
    ce_dataset = xr.Dataset(
        {
            "CE": (
                ("j", "i"),
                grid,
            )
        },
        attrs={
            "units": "-",
            "long_name": "Whole squares",
            "name": "CE",
        },
    )
    return ds.assign(CE=ce_dataset["CE"])
