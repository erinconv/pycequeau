"""
This is a helper to download and manage all the copernicus files
"""

import os
import re
import sys
from typing import Any, Union
import subprocess
from urllib.parse import quote
import xml.etree.ElementTree as ET
import requests
import xarray as xr
import numpy as np
import rioxarray as rxr
import xrspatial as xrs
from osgeo import gdal

from .manage_files import delete_raster_and_sidecars

class CopernicusDEMProcessor:
    """_summary_
    """
    def __init__(self, project_path: str, extent: Union[list, tuple]):
        self._project_path = project_path
        self._project_extent = extent
        self._cop_files_path = os.path.join(self._project_path, "COP_files")
        self._merged_dem = os.path.join(self._project_path, "geographic/merged_DEM.tif")
        self._merged_wbm = os.path.join(self._project_path, "geographic/merged_WBM.tif")

        # Clipped files lists
        self._clipped_merged_dem = os.path.join(self._project_path, "geographic/DEM.tif")
        self._clipped_merged_wbm = os.path.join(self._project_path, "geographic/WBM.tif")

    def download_and_merge(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        self._check_project_structure()
        if os.path.exists(self._merged_wbm) and os.path.exists(self._merged_dem):
            print("COP dem and wbm files already exist for this project")
        else:
            self._download_COP_files()
            self._merge_rasters()
            self._clip_rasters()
        self._clean_files()
        return 0

    def _clean_files(self):
        """Remove merged raster files and their GDAL sidecars after clipping."""
        files_to_clean = [self._merged_dem, self._merged_wbm]
        for path in files_to_clean:
            delete_raster_and_sidecars(path)

    def _check_project_structure(self):
        """_summary_

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        local_dirs = os.listdir(self._project_path)
        if not "geographic" in local_dirs:
            os.mkdir(os.path.join(self._project_path, "geographic"))
        if not "meteo" in local_dirs:
            os.mkdir(os.path.join(self._project_path, "meteo"))
        if not "results" in local_dirs:
            os.mkdir(os.path.join(self._project_path, "results"))

    # Setter to define the extent of the project
    @property
    def project_extent(self):
        """project_extent getter
        """
        return self._project_extent

    @project_extent.setter
    def project_extent(self, extent: tuple):
        self._project_extent = extent

    def _list_s3_prefixes(self, bucket: str, prefix: str, delimiter: str = '/') -> list:
        """List all prefixes (folders) in an S3 bucket using REST API."""
        base_url = f"https://{bucket}.s3.eu-central-1.amazonaws.com/"
        # S3 REST API parameters
        params = {
            'list-type': '2',
            'prefix': prefix,
            'delimiter': delimiter,
            'max-keys': '1000'}
        prefixes = []
        continuation_token = None
        while True:
            if continuation_token:
                params['continuation-token'] = continuation_token
            # Build query string
            query_string = '&'.join([f"{k}={quote(v)}" for k, v in params.items() if v])
            url = f"{base_url}?{query_string}"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            # Parse XML response (S3 returns XML)
            root = ET.fromstring(response.text)
            # Find all CommonPrefixes
            ns = {'': 'http://s3.amazonaws.com/doc/2006-03-01/'}
            for prefix_elem in root.findall('.//CommonPrefixes/Prefix', ns):
                prefixes.append(prefix_elem.text)
            # Check for continuation token
            is_truncated = root.find('.//IsTruncated', ns)
            if is_truncated is None or is_truncated.text != 'true':
                break
            next_token = root.find('.//NextContinuationToken', ns)
            if next_token is not None:
                continuation_token = next_token.text
        return prefixes

    def _parse_tile_coordinates(self, prefix: str):
        """
        Parse tile coordinates from Copernicus DSM folder prefix.
        Format: Copernicus_DSM_COG_10_N{lat}_{dec}_E{lon}_{dec}_DEM/
        Example: Copernicus_DSM_COG_10_N00_00_E006_00_DEM/ -> N00.00, E006.00
        Or: Copernicus_DSM_COG_10_N48_00_W070_00_DEM/ -> N48.00, W070.00
        """
        # Remove trailing slash
        base = prefix.rstrip('/')
        # Try to find coordinate patterns
        # Look for patterns like N48_00_E006_00 or N00_00_E006_00
        lon_pattern = r'([EW])(\d{3})_(\d{2})'
        lat_pattern = r'([NS])(\d{2})_(\d{2})'
        lon_match = re.search(lon_pattern, base)
        lat_match = re.search(lat_pattern, base)
        if lon_match and lat_match:
            lon_dir = lon_match.group(1)
            lon_deg = int(lon_match.group(2))
            lon_dec = int(lon_match.group(3))
            lat_dir = lat_match.group(1)
            lat_deg = int(lat_match.group(2))
            lat_dec = int(lat_match.group(3))
            # Convert to decimal degrees
            lon = lon_deg + lon_dec / 100.0
            lat = lat_deg + lat_dec / 100.0
            # Apply direction
            if lon_dir == 'W':
                lon = -lon
            if lat_dir == 'S':
                lat = -lat
            return lon, lat
        else:
            return None, None

    def _tile_intersects_bounds(self, 
                               tile_lon: float,
                               tile_lat: float,
                               min_lon: float,
                               max_lon: float,
                               min_lat: float,
                               max_lat: float):
        """
        Check if a tile intersects with the bounding box.
        Assume tiles are 1x1 degree.
        """
        # Tile covers [tile_lon, tile_lon+1] x [tile_lat, tile_lat+1]
        tile_max_lon = tile_lon + 1
        tile_max_lat = tile_lat + 1
        # Check for intersection
        return not (tile_max_lon < min_lon or tile_lon > max_lon or 
                    tile_max_lat < min_lat or tile_lat > max_lat)

    def _list_s3_objects(self, bucket: str, prefix: str):
        """List all objects in an S3 bucket prefix using REST API."""
        base_url = f"https://{bucket}.s3.eu-central-1.amazonaws.com/"
        
        params = {
            'list-type': '2',
            'prefix': prefix,
            'max-keys': '1000'
        }
        
        objects = []
        try:
            # Build query string
            query_string = '&'.join([f"{k}={quote(v)}" for k, v in params.items()])
            url = f"{base_url}?{query_string}"
            response = requests.get(url, timeout=50)
            response.raise_for_status()
            # Parse XML response
            root = ET.fromstring(response.text)
            ns = {'': 'http://s3.amazonaws.com/doc/2006-03-01/'}
            for content in root.findall('.//Contents', ns):
                key = content.find('Key', ns).text
                size = int(content.find('Size', ns).text)
                objects.append((key, size))
        
        except Exception as e:
            print(f"Error listing objects: {e}")
            return []
        return objects

    def _download_s3_object(self, bucket: str, key: str, output_path: str):
        """Download an object from S3 using REST API."""
        url = f"https://{bucket}.s3.eu-central-1.amazonaws.com/{key}"
        try:
            response = requests.get(url, stream=True, timeout=300)
            response.raise_for_status()
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            return True
        except Exception as e:
            print(f"Error downloading: {e}")
            return False

    def _download_COP_files(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        min_lon = self._project_extent[0]
        max_lon = self._project_extent[1]
        min_lat = self._project_extent[2]
        max_lat = self._project_extent[3]
        matched_tiles = []
        all_tiles_info = []  # Store all tiles for debugging
        prefixes = self._list_s3_prefixes('copernicus-dem-30m', 'Copernicus_DSM_COG_10')
        for prefix in prefixes:
            tile_lon, tile_lat = self._parse_tile_coordinates(prefix)
            all_tiles_info.append((prefix, tile_lon, tile_lat))
            if tile_lon is not None and tile_lat is not None:
                if self._tile_intersects_bounds(tile_lon, tile_lat, min_lon, max_lon, min_lat, max_lat):
                    matched_tiles.append((prefix, tile_lon, tile_lat))
        print(f"\n{'='*70}")
        print(f"Matched {len(matched_tiles)} tiles for your region:")
        print(f"{'='*70}")
        # Make sure that the COP_DEM folder exist inside the geographic folder
        cop_files = os.path.join(self._project_path, "geographic", "COP_files")
        os.makedirs(cop_files, exist_ok=True)
        if matched_tiles:
            for prefix, tile_lon, tile_lat in matched_tiles:
                # Get the main DEM file and water body mask from this tile folder
                # The main file is: prefix/Copernicus_DSM_COG_10_..._DEM.tif
                # Water body mask is: prefix/Copernicus_DSM_COG_10_..._WBM.tif
                # List files in this folder to find the main DEM and water body mask
                tile_objects = self._list_s3_objects('copernicus-dem-30m', prefix)
                if tile_objects:
                    for key, size in tile_objects:
                        # Skip aux files and previews
                        if 'PREVIEW' in key:
                            continue
                        
                        filename = key.split('/')[-1]
                        
                        # Download DEM file
                        if key.endswith('_DEM.tif'):
                            output_path = os.path.join(cop_files, filename)
                            # Skip if file already exists
                            if os.path.exists(output_path):
                                print(f"OK {filename} (lon: {tile_lon:.2f}, lat: {tile_lat:.2f}) - already exists, skipping")
                                continue
                            print(f"Downloading DEM: {filename} (lon: {tile_lon:.2f}, lat: {tile_lat:.2f}, size: {size:,} bytes)...")
                            if self._download_s3_object('copernicus-dem-30m', key, output_path):
                                print(f"  Successfully downloaded to: {output_path}")
                            else:
                                print("  Error downloading DEM file")
                        
                        # Download water body mask file (common patterns: _WBM.tif, _WATER.tif, _WATER_BODY.tif)
                        elif key.endswith('_WBM.tif') or key.endswith('_WATER.tif') or key.endswith('_WATER_BODY.tif'):
                            output_path = os.path.join(cop_files, filename)
                            # Skip if file already exists
                            if os.path.exists(output_path):
                                print(f"OK WBM {filename} (lon: {tile_lon:.2f}, lat: {tile_lat:.2f}) - already exists, skipping")
                                continue
                            print(f"Downloading Water Body Mask: {filename} (lon: {tile_lon:.2f}, lat: {tile_lat:.2f}, size: {size:,} bytes)...")
                            if self._download_s3_object('copernicus-dem-30m', key, output_path):
                                print(f"  Successfully downloaded to: {output_path}")
                            else:
                                print("  Error downloading water body mask file")
        return 0

    def _merge_rasters(self):
        """Merge downloaded DEM and WBM raster tiles into single files using gdal_merge."""
        # Find gdal_merge.py - try multiple possible locations
        possible_paths = [
            os.path.join(os.path.dirname(gdal.__file__), 'scripts', 'gdal_merge.py'),
            os.path.join(os.path.dirname(os.path.dirname(gdal.__file__)), 'osgeo_utils', 'gdal_merge.py'),
            'gdal_merge.py'  # If in PATH
        ]

        gdal_scripts_path = None
        for path in possible_paths:
            if os.path.exists(path) or path == 'gdal_merge.py':
                gdal_scripts_path = path
                break

        if gdal_scripts_path is None:
            print("Error: Could not find gdal_merge.py script")
            return -1

        # List all files stored in the COP_files folders
        cop_files = os.path.join(self._project_path, "geographic", "COP_files")
        output_path = os.path.join(self._project_path, "geographic")

        if not os.path.exists(cop_files):
            print(f"Error: Directory {cop_files} does not exist")
            return -1

        files_list = os.listdir(cop_files)
        # Merged files in a list
        # self._merged_files = []
        # Process each data type separately (DEM and WBM)
        for data_type in ["DEM", "WBM"]:
            raster_list = []
            for file in files_list:
                if file.endswith(f"{data_type}.tif"):
                    raster_list.append(os.path.join(cop_files, file))

            if not raster_list:
                print(f"No {data_type} files found to merge")
                continue

            # Create output filename
            output_file = os.path.join(output_path, f"merged_{data_type}.tif")

            # Remove the output file if it exists
            if os.path.exists(output_file):
                os.remove(output_file)

            # Build command: python gdal_merge.py -o output_file input1.tif input2.tif ...
            # Or: gdal_merge -o output_file input1.tif input2.tif ... (if in PATH)
            if gdal_scripts_path == 'gdal_merge.py':
                cmd = ['gdal_merge', '-o', output_file] + raster_list
            else:
                cmd = [sys.executable, gdal_scripts_path, '-o', output_file] + raster_list

            try:
                print(f"Merging {len(raster_list)} {data_type} files into {output_file}...")
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"Successfully merged {data_type} files")
            except subprocess.CalledProcessError as e:
                print(f"Error merging {data_type} files: {e.stderr}")
                return -1

        return 0

    def _clip_rasters(self):
        """Clip merged DEM and WBM rasters to project extent using gdal raster clip."""
        if self._project_extent is None:
            print("Error: Project extent not set. Cannot clip rasters.")
            return -1

        # Extract bounding box from project extent
        # Format: (min_lon, max_lon, min_lat, max_lat)
        min_lon, max_lon, min_lat, max_lat = self._project_extent

        # GDAL raster clip expects bbox as: xmin,ymin,xmax,ymax (lon,lat order)
        bbox_str = f"{min_lon},{min_lat},{max_lon},{max_lat}"

        geographic_path = os.path.join(self._project_path, "geographic")

        # Process each data type separately (DEM and WBM)
        for data_type in ["DEM", "WBM"]:
            input_file = os.path.join(geographic_path, f"merged_{data_type}.tif")
            output_file = os.path.join(geographic_path, f"{data_type}.tif")

            # Check if merged file exists
            if not os.path.exists(input_file):
                print(f"Warning: Merged {data_type} file not found: {input_file}")
                continue

            # Build command: gdal raster clip --bbox=xmin,ymin,xmax,ymax --bbox-crs=EPSG:4326 input.tif output.tif --overwrite
            cmd = [
                'gdal', 'raster', 'clip',
                '--bbox', bbox_str,
                '--bbox-crs', 'EPSG:4326',  # Assuming extent is in WGS84 lat/lon
                input_file,
                output_file,
                '--overwrite'
            ]

            try:
                print(f"Clipping {data_type} raster to project extent ({bbox_str})...")
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"Successfully clipped {data_type} to: {output_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error clipping {data_type}: {e.stderr}")
                return -1
            except FileNotFoundError:
                print("Error: 'gdal' command not found. Make sure GDAL is installed and in PATH.")
                return -1

        return 0

    @classmethod
    def basic_dem_condition(cls, in_dem: str, in_wbm: str):
        """_summary_

        Returns:
            _type_: _description_
        """
        base_dir = os.path.dirname(in_dem)
        # Open rasters
        dem_raster = rxr.open_rasterio(in_dem, chunks="auto").squeeze()
        wbm_raster = rxr.open_rasterio(in_wbm, chunks="auto").squeeze()

        # Get ocean pixels if they exist
        ocean_pixels = wbm_raster.where(wbm_raster == 1)
        # Set ocean pixels to no data in the DEM
        dem_raster = dem_raster.where(ocean_pixels.isnull(), np.nan)

        # Now compute the surface water body areas as boolean mask
        boolean_surface_water = xr.where(wbm_raster > 1,1,np.nan)
        boolean_surface_water = xr.where(np.isnan(boolean_surface_water),1,0)

        # Compute distance allocation on the surface water body areas
        distance = xrs.proximity(boolean_surface_water, distance_metric='GREAT_CIRCLE')

        # Reduce the proximity values by logartitmic transformation
        distance = np.log(distance + 1)/10 - 2 #Add negative 2 to make sure water areas are lower than the surroundings.

        # Substract the distance from the DEM
        no_flat_dem = dem_raster - distance
        no_flat_dem_path = os.path.join(base_dir, "DEM_conditioned.tif")

        # distance = boolean_surface_water
        # save dataset as tif
        dem_crs = dem_raster.rio.crs
        dem_dtype = dem_raster.dtype
        # Replace the no data value in the DEM to -99999
        dem_nd = -99999
        dem_raster = dem_raster.where(np.isnan(dem_raster), dem_nd)
        no_flat_dem.rio.write_crs(dem_crs, inplace=True)
        no_flat_dem.rio.write_nodata(dem_nd, inplace=True)
        no_flat_dem.rio.to_raster(no_flat_dem_path, dtype=dem_dtype)
        # save dataset as xarray
        return 0