from __future__ import annotations

import os
import json
import numpy as np
import geopandas as gpd
import xarray as xr
import pandas as pd
from shapely.geometry import Point, LineString
from src.physiographic.base import Basin
from src.core import projections
# import shapely.geometry


class Parameters:
    """_summary_
    """

    def __init__(self, bassinVersant: Basin) -> None:
        """_summary_

        Args:
            bassinVersant (Basin): _description_
        """
        self.ctp = {0}
        self.lac = 0
        self.surface = 0
        self.basin_structure = bassinVersant
        self.parametres = None
        self.option = None
        self.sol = None
        self.solInitial = None
        self.transfert = None
        self.fonte = None
        self.evapo = None
        self.dli = None
        self.qualite = None
        self.main_channel_length = None
        self.time_of_concentrations = None
        self.slope_channel = None
        self._jonei = None

    def day_max_insolation(self, meteo_file_name: str):
        """_summary_

        Args:
            meteo_file_name (str): _description_

        Returns:
            _type_: _description_
        """
        meteo_file_name = os.path.join(self.basin_structure.project_path,
                                       "meteo", meteo_file_name)
        # Open the meteo file name
        nc_meteo = xr.open_dataset(meteo_file_name)
        # Get the variable tmax in the first cp
        tmax = nc_meteo["tMax"].to_pandas()
        # tmax = nc_meteo["rayonnement"][:,0].values
        time = pd.to_datetime(nc_meteo["pasTemp"].values-719529, unit='D')
        tmax.index = time
        # Get the julian days
        multi_annual = tmax.groupby(tmax.index.day_of_year).mean()
        idx_max_temp = multi_annual.idxmax()
        jour = multi_annual.index[idx_max_temp].values
        jonei = np.mean(jour - 365/4)
        self._jonei = int(jonei)

    def create_parameter_structure(self):
        """_summary_
        """
        self.set_longwave_radiation_parameters()
        self.parametres = {"option": self.option,
                           "sol": self.sol,
                           "solInitial": self.solInitial,
                           "transfert": self.transfert,
                           "ctp": self.ctp,
                           "lac": self.lac,
                           "surface": self.surface,
                           "fonte": self.fonte,
                           "evapo": self.evapo,
                           "qualite": self.qualite,
                           "dli": self.dli}
        with open(os.path.join(self.basin_structure.project_path,
                               "results", "parameters.json"), "w",
                  encoding="utf-8") as outfile:
            json.dump(self.parametres, outfile, indent=4, default=tuple)
        outfile.close()

    def set_option(self, values: np.ndarray):
        """_summary_

        Args:
            values (np.ndarray): _description_
        """
        # The default values can be changed using the method set_maximum_insolation_day
        self.option = {"ipassim": 24,
                       "moduleFonte": values[0],
                       "moduleEvapo": values[1],
                       "moduleOmbrage": 0,
                       "moduleDLI": 1,
                       "calculQualite": values[2],
                       "jonei": self._jonei,
                       "joeva": self._jonei}

    def set_sol(self, values: np.ndarray):
        """_summary_

        Args:
            values (np.ndarray): _description_
        """
        # Open the carreaux entier
        carreux_entier_name = os.path.join(self.basin_structure.project_path,
                                           "results",
                                           "carreauxEntiers.csv")
        CE_df = pd.read_csv(carreux_entier_name, index_col=0)
        impermeable = CE_df["pctImpermeable"].to_list()
        self.sol = {"cin_s": values[0],
                    "cvmar": values[1],
                    "cvnb_s": values[2],
                    "cvnh_s": values[3],
                    "cvsb": values[4],
                    "cvsi_s": values[5],
                    "xinfma": values[6],
                    "hinf_s": values[7],
                    "hint_s": values[8],
                    "hmar": values[9],
                    "hnap_s": values[10],
                    "hpot_s": values[11],
                    "hsol_s": values[12],
                    "hrimp_s": values[13],
                    "tri_s": 0,
                    "xla": self._compute_xla(),
                    }

    def _compute_xla(self) -> int:
        """
        The format for the mean latitudde for the CEQUEAU model is given as follows:
        let us take a latitude float coordinate:
        lat = XX.xx
        where XX is the whole part and xx are the two first decimal positions. The xla parameter
        inside the CEQUEAU mudel must be given as an integer number in the next format:
        xla = XXxx
        """
        # Open the watershed shp file as geopandas
        watershed = gpd.read_file(
            self.basin_structure.Basin.replace(".tif", ".shp"))
        centroid = watershed.centroid
        # Get the epsg code from the shp
        epsg_code = centroid.crs.srs
        x = centroid.x.values.tolist()
        y = centroid.y.values.tolist()
        # Convert the centroid from utm to lat lon
        x, y = projections.utm_to_latlon(x, y, epsg_code)
        # Convert the y to a string to a string and round irt to two decimals
        y_str = str(round(y[0], 2))
        y_str = y_str.replace(".", "")
        if len(y_str) < 4:
            y_str = y_str + "0"
        # Convert it back to an integer number
        y_int = int(y_str)
        return y_int

    def set_solinitial(self, values: np.ndarray):
        """_summary_

        Args:
            values (np.ndarray): _description_
        """
        self.solInitial = {"hsini": values[0],
                           "hnini": values[1],
                           "hmini": values[2],
                           "q0": values[3],
                           "tmur": values[4],
                           "tstock": values[5],
                           }

    def set_transfert(self, values: np.ndarray):
        """_summary_

        Args:
            values (np.ndarray): _description_
        """

        # self.create_cequeau_stream_network()
        self.transfert = {"exxkt": values[0],
                          "zn": self.compute_tc(),
                          "tc_struct": self.time_of_concentrations}

    def compute_tc(self):
        r"""
        All the methods used in this program are based in the equations presented in:
        https://www.tandfonline.com/doi/pdf/10.1080/02626667.2011.644244

        and here:
        https://www.redalyc.org/articulo.oa?id=49622372006
        """
        cp_struct_name = os.path.join(self.basin_structure.project_path,
                                      "results",
                                      "carreauxPartiels.csv")
        carreaux_partiels = pd.read_csv(cp_struct_name, index_col=0)
        # find the difference between the highest and lowest point in altitute
        hmin = carreaux_partiels["altitudeMoy"].min()
        hmax = carreaux_partiels["altitudeMoy"].max()
        hmean = carreaux_partiels["altitudeMoy"].mean()
        # Get the basin area in km2
        basin_area = self.basin_structure.bassinVersant.get(
            "superficieCE")*carreaux_partiels["pctSurface"].sum()/100
        # Create time of concentration structure
        self.time_of_concentrations = {"Kirpich": self.tc_kirpich(self.main_channel_length, hmax - hmin),
                                       "Giandotti": self.tc_giandotti(basin_area, self.main_channel_length, hmean - hmin),
                                       "Department_of_public_work": self.tc_pw(self.main_channel_length, hmax - hmin)}
        # Convert to df
        df = pd.DataFrame.from_dict(self.time_of_concentrations,
                                    orient="index")
        # Compute average Tc
        avg_tc = df[0].mean(axis=0)

        return avg_tc

    # @classmethod
    # def tc_clark(cls, basin_area: float,
    #              main_channel_length: float,
    #              height_differences: float) -> float:
    #     """
    #     This module uses the Department of Public Works (1995) formulation to compute the time of concentration:

    #     Raises:
    #         ValueError: _description_

    #     Returns:
    #         _type_: _description_
    #     """
    #     slope_channel = height_differences/main_channel_length
    #     tc = 0.335*(basin_area/np.sqrt(slope_channel))**0.593
    #     tc = tc/24
    #     return 0

    @classmethod
    def tc_pw(cls, main_channel_length: float,
              height_differences: float):
        r"""
        This module uses the Department of Public Works (1995) formulation to compute the time of concentration:

        .. math::

            T_{c} = 60 \left ( 11.9 \frac{L^{3}}{H} \right )^{0.385}

        where:
            - :math:`L` is the length of the main channel in :math:`mi`
            - :math:`H` is the maximun elevation difference in :math:`ft`
            - :math:`T_{c}` is the time of concentration in :math:`min`

        Args:
            main_channel_length (float): _description_
            height_differences (float): _description_

        Returns:
            _type_: _description_
        """
        tc = 60*(11.9*(0.000621371*main_channel_length)
                 ** 3/(3.28084*height_differences))**0.385
        # Convert to days
        tc = tc/1440
        return tc

    @classmethod
    def tc_giandotti(cls, basin_area: float,
                     main_channel_length: float,
                     height_differences: float) -> float:
        r"""
        This module uses the Giandotti (1934) formulation to compute the time of concentration:

        .. math::

            T_{c} = \frac{4\sqrt{A} + 1.5L}{0.8\sqrt{H}}

        where:
            - :math:`A` is the area of the basin in :math:`\text{km}^{2}`
            - :math:`L` is the length of the main channel in :math:`\text{km}`
            - and :math:`H` is the the difference between the mean basin elevation and the outlet elevation
            - :math:`T_{c}` is the time of concentration in hours.

        This formula was calibrated on 12 basins with drainage areas between 170 and 70 000 :math:`km^{2}`

        Args:
            basin_area (float): _description_
            main_channel_length (float): _description_
            height_differences (float): _description_

        Returns:
            _type_: _description_
        """
        # Compute tc in hours
        tc = (4*np.sqrt(basin_area)+1.5*main_channel_length/1e3) / \
            (0.8*np.sqrt(height_differences))
        # Convert tc to days
        tc = tc/24
        return tc

    @classmethod
    def tc_kirpich(cls, main_channel_length: float,
                   height_differences: float) -> float:
        r"""
        This module uses the Kirpich (1940) formulation to compute the time of concentration:

        .. math::

            T_{c} = 0.0078 L^{0.77}S^{-0.385}

        Where :math:`L` is the length of the main channel (ft) in :math:`m`, :math:`S` is the mean
        slope of the basin in  :math:`\text{m m}^{-1}` which is computed as follows:

        .. math::

            S = \frac{H_{max} - H_{min}}{L}

        :math:`T_{c}` is given in minutes which is subsequently converted into days

        Args:
            main_channel_length (float): _description_
            carreaux_partiels (pd.DataFrame): _description_

        Returns:
            _type_: _description_
        """
        s_basin = height_differences/main_channel_length
        # The time of concentration is given in minutes
        tc = 0.0078*(3.28084*main_channel_length)**(0.77)*(s_basin**(-0.385))
        # Convert to days
        tc = tc/1440
        return tc

    def create_cequeau_stream_network(self, area_th=0.01) -> float:
        """_summary_

        Args:
            area_th (float, optional): _description_. Defaults to 0.01.

        Returns:
            float: _description_
        """
        # Open the
        CP_fishnet = gpd.read_file(self.basin_structure.cp_fishnet_name)
        cp_struct_name = os.path.join(self.basin_structure.project_path,
                                      "results",
                                      "carreauxPartiels.csv")
        carreaux_partiels = pd.read_csv(cp_struct_name, index_col=0)
        CP_fishnet["x_c"] = CP_fishnet.centroid.x
        CP_fishnet["y_c"] = CP_fishnet.centroid.y

        # Array to store the line shape
        lines = np.zeros(len(CP_fishnet), dtype=LineString)
        # Mask small streams
        max_area = CP_fishnet["cumulArea"].max()
        idx_small_area = CP_fishnet["cumulArea"] > area_th*max_area
        CP_fishnet.index = CP_fishnet["newCPid"].values
        carreaux_partiels.index = CP_fishnet["newCPid"].values
        df_line = pd.DataFrame(columns=["names", "cumulArea", "slope"],
                               data=np.empty(shape=[len(CP_fishnet), 3]),
                               index=CP_fishnet.index)
        df_line["cumulArea"] = CP_fishnet["cumulArea"].values
        for i, _ in CP_fishnet.iterrows():
            cp_aval = carreaux_partiels.loc[i, "idCPAval"]
            if cp_aval == 0:
                p1 = Point(CP_fishnet.loc[i, "x_c"],
                           CP_fishnet.loc[i, "y_c"])
                p2 = Point(CP_fishnet.loc[i, "x_c"],
                           CP_fishnet.loc[i, "y_c"])
                lines[i-1] = LineString([p2, p1])
                # carreaux_partiels["altitudeMoy"]
                df_line.loc[i, "names"] = f"CP{i} to CP{cp_aval}"
                slope = 0
                df_line.loc[i, "slope"] = slope
            else:
                p1 = Point(CP_fishnet.loc[i, "x_c"],
                           CP_fishnet.loc[i, "y_c"])
                p2 = Point(CP_fishnet.loc[cp_aval, "x_c"],
                           CP_fishnet.loc[cp_aval, "y_c"])
                lines[i-1] = LineString([p2, p1])
                df_line.loc[i, "names"] = f"CP{i} to CP{cp_aval}"
                dh = carreaux_partiels.loc[i, "altitudeMoy"] - \
                    carreaux_partiels.loc[cp_aval, "altitudeMoy"]
                df_line.loc[i, "slope"] = abs(dh)/lines[i-1].length

        # df_line = df_line.loc[idx_small_area.values]
        gpdf_line = gpd.GeoDataFrame(df_line,
                                     geometry=lines,
                                     crs=CP_fishnet.geometry.crs)
        gpdf_line = gpdf_line.loc[idx_small_area.values]
        # Open the outlet routes
        outlet_file = os.path.join(self.basin_structure.project_path,
                                   "results",
                                   "outlet_routes.csv")
        outlet_routes = pd.read_csv(outlet_file, header=None)
        # find all the complete routes
        idy, = np.where(outlet_routes.iloc[:, -1] != 0)
        outlet_routes_complete = outlet_routes.iloc[idy, :]
        # Compute the length of all the routes
        lengths = []
        # Find the largest CP with stream
        idx_max = gpdf_line.index.max()
        for i in range(len(outlet_routes_complete)):
            routes = np.array(outlet_routes_complete.iloc[i, :])
            routes = routes[routes < idx_max]
            sub_gpdf = gpdf_line.loc[routes, :]
            lengths.append(sub_gpdf.geometry.length.sum())
        # Find the maximun length value
        idx_longest_path = lengths.index(max(lengths))
        main_route = np.array(outlet_routes_complete.iloc[idx_longest_path, :])
        main_route = main_route[main_route < idx_max]
        gpdf_line["main_path"] = 0
        gpdf_line.loc[main_route, "main_path"] = 1
        # Open the outlet routes
        streams_file = os.path.join(self.basin_structure.project_path,
                                    "geographic",
                                    "streams_cequeau.shp")
        gpdf_line.to_file(streams_file)
        self.main_channel_length = max(lengths)
        self.slope_channel = gpdf_line["slope"].where(
            gpdf_line["slope"] > 0).mean()

    def set_fonte(self, values: np.ndarray, model: int):
        # Parameters for the cequeau model
        self.fonte = {}
        self.fonte["cequeau"] = {}
        # DEGREE-DAY
        self.fonte["cequeau"]["strne_s"] = values[0]
        self.fonte["cequeau"]["tfc_s"] = values[1]
        self.fonte["cequeau"]["tfd_s"] = values[2]
        self.fonte["cequeau"]["tsc_s"] = values[3]
        self.fonte["cequeau"]["tsd_s"] = values[4]
        self.fonte["cequeau"]["ttd"] = values[5]
        self.fonte["cequeau"]["tts_s"] = values[6]
        self.fonte["cequeau"]["jonei"] = self.option["jonei"]
        self.fonte["cequeau"]["tmur"] = self.solInitial["tmur"]
        self.fonte["cequeau"]["tstock"] = self.solInitial["tstock"]

        # UEB
        self.fonte["UEB"] = {}
        self.fonte["UEB"]["strne_s"] = values[0]
        self.fonte["UEB"]["K_s"] = 0.15
        self.fonte["UEB"]["z0"] = 0.003
        self.fonte["UEB"]["aep"] = 0.2
        self.fonte["UEB"]["K_sat"] = 350
        self.fonte["UEB"]["rho_s"] = 450
        self.fonte["UEB"]["melt_frac"] = 0.99
        self.fonte["UEB"]["melt_thr"] = 0
        self.fonte["UEB"]["hours"] = 16
        self.fonte["UEB"]["z"] = 2.0
        self.fonte["UEB"]["avo"] = 0.8
        self.fonte["UEB"]["airo"] = 0.6
        self.fonte["UEB"]["Lc"] = 0.05
        self.fonte["UEB"]["fstab"] = 1
        self.fonte["UEB"]["D"] = 0.001
        self.fonte["UEB"]["de"] = 0.4
        self.fonte["UEB"]["snow_temp_method"] = 1
        # Initial conditions
        self.fonte["UEB"]["w"] = 0
        self.fonte["UEB"]["ub"] = -2000
        self.fonte["UEB"]["E"] = 0
        self.fonte["UEB"]["tausn"] = 0
        self.fonte["UEB"]["tsurf"] = 0
        self.fonte["UEB"]["tave"] = 0
        self.fonte["UEB"]["Mr"] = 0
        self.fonte["UEB"]["albedo"] = 0.25
        # CEMANAIEGE
        self.fonte["cemaNeige"] = {}
        self.fonte["cemaNeige"]["strne"] = values[0]
        self.fonte["cemaNeige"]["Kf"] = 15
        self.fonte["cemaNeige"]["Tf"] = 0
        self.fonte["cemaNeige"]["CTg"] = 0.85
        self.fonte["cemaNeige"]["theta"] = 0.8
        self.fonte["cemaNeige"]["Gseuil"] = 250*0.9
        self.fonte["cemaNeige"]["Vmin"] = 0.5
        self.fonte["cemaNeige"]["Zmed"] = 300
        # Initial conditions
        self.fonte["cemaNeige"]["eTg"] = 0.1
        self.fonte["cemaNeige"]["G"] = 0

        a = 1

    def set_evapo(self, values: np.ndarray, model: int):
        # CEQUEAU - THORNWAIT
        if model == 1:
            self.evapo = {"cequeau": {
                "joeva": self.option["joeva"],
                "evnap": values[0],
                "xaa": values[1],
                "xit": values[2]
            }
            }
        # KPENNMAN
        elif model == 2:
            self.evapo = {"evnap": values[0]}
        # PriestleyTaylor
        elif model == 3:
            self.evapo = {"alpha": values[0],
                          "evnap": values[1]}
        # McGuinness
        elif model == 4:
            self.evapo = {"evnap": values[0]}
        # PenmanMont
        elif model == 5:
            self.evapo = {"evnap": values[0]}
        # Morton
        elif model == 6:
            self.evapo = {"alpha": values[0],
                          "evnap": values[1]}

    def set_longwave_radiation_parameters(self):
        self.dli = {}
        self.dli["m1"] = {}
        self.dli["m2"] = {}
        self.dli["m3"] = {}
        self.dli["m4"] = {}
        self.dli["m5"] = {}
        self.dli["m6"] = {}
        self.dli["m7"] = {}
        self.dli["m8"] = {}

        # Fill the parameters with default values
        self.dli["m1"] = {}
        self.dli["m1"]["u"] = 0
        self.dli["m1"]["v"] = 1
        self.dli["m1"]["a"] = 0.8171
        # m2
        self.dli["m2"]["u"] = 0.2184
        self.dli["m2"]["v"] = 2.4311
        self.dli["m2"]["a"] = 9.3645e-6
        # m3
        self.dli["m3"]["u"] = 0.17
        self.dli["m3"]["v"] = 4
        self.dli["m3"]["a"] = 0.2610
        self.dli["m3"]["b"] = -7.77e-2
        # m4
        self.dli["m4"]["u"] = 0.17
        self.dli["m4"]["v"] = 2
        self.dli["m4"]["a"] = 0.74
        self.dli["m4"]["b"] = 0.0065
        # m5
        self.dli["m5"]["u"] = 0.2187
        self.dli["m5"]["v"] = 1.6689
        self.dli["m5"]["a"] = 1.24
        self.dli["m5"]["b"] = 0.1429
        # m6
        self.dli["m6"]["u"] = 0.1296
        self.dli["m6"]["v"] = 4
        self.dli["m6"]["a"] = 1.08
        self.dli["m6"]["b"] = 2016
        # m7
        self.dli["m7"]["u"] = 0.1827
        self.dli["m7"]["v"] = 3.4472
        self.dli["m7"]["a"] = 46.5
        self.dli["m7"]["b"] = 1.2
        self.dli["m7"]["c"] = 3
        self.dli["m7"]["d"] = 0.5
        # m8
        self.dli["m8"]["u"] = 0.1533
        self.dli["m8"]["v"] = 4
        self.dli["m8"]["a"] = 0.72
        self.dli["m8"]["b"] = 0.0090
        self.dli["m8"]["c"] = 0.078
        return 0

    def set_qualite(self, values: np.ndarray, model=1):
        # CEQUEAU model
        if model == 1:
            temperat = {"crayso": values[2],
                        "crayin": values[3],
                        "cevapo": values[4],
                        "cconve": values[5],
                        "crigel": values[6],
                        "tnap": values[7],
                        "bassol": values[8],
                        "corsol": values[9],
                        # This is true if the simulations starts in the winter
                        "panap": 1,
                        # This is true if the simulations starts in the winter
                        "tinit": 0}
            self.qualite = {"cequeau": {
                "coprom": values[0],
                            "colarg": values[1],
                            "temperat": temperat
                            }
                            }
        # TODO: In the future, new water temperature models can be added
        # Future model?
        elif model != 1:
            raise ValueError(
                "No model with this label has ben yet added to the CEQUEAU model")

    # TODO: Add the parameters for the barrage
    def set_barrage(self, values: np.ndarray):
        pass

# Meleze_Struct.parametres.option,
# Meleze_Struct.parametres.sol,
# Meleze_Struct.parametres.solInitial, Meleze_Struct.parametres.transfert, Meleze_Struct.parametres.ctp, Meleze_Struct.parametres.lac, Meleze_Struct.parametres.surface, Meleze_Struct.parametres.interpolation, Meleze_Struct.parametres.fonte, Meleze_Struct.parametres.evapo, Meleze_Struct.parametres.neige, Meleze_Struct.parametres.qualite, Meleze_Struct.parametres.barrage
