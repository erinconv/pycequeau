from __future__ import annotations

import numpy as np
import os
import json
import geopandas as gpd
import pandas as pd
import shapely.geometry
from shapely.geometry import Point, LineString
from src.physiographic.base import Basin
from src.core import projections


class Parameters:
    def __init__(self, bassinVersant: Basin) -> None:
        self.ctp = {0}
        self.lac = 0
        self.surface = 0
        self.basin_structure = bassinVersant
        pass

    def create_parameter_structure(self):
        self.parametres = {"option": self.option,
                           "sol": self.sol,
                           "solInitial": self.solInitial,
                           "transfert": self.transfert,
                           "ctp": self.ctp,
                           "lac": self.lac,
                           "surface": self.surface,
                           "fonte": self.fonte,
                           "evapo": self.evapo,
                           "qualite": self.qualite}
        with open(os.path.join(self.basin_structure._project_path, "results", "parameters.json"), "w") as outfile:
            json.dump(self.parametres, outfile, indent=4, default=tuple)
        outfile.close()

    def set_option(self, values: np.ndarray):
        # The default values can be changed using the method set_maximum_insolation_day
        self.option = {"ipassim": 24,
                       "moduleFonte": values[0],
                       "moduleEvapo": values[1],
                       "calculQualite": values[2],
                       "jonei": 80,
                       "joeva": 80}

    def set_maximum_insolation_day(self, jonei: int):
        self.option["jonei"] = jonei
        self.option["joeva"] = jonei

    def set_sol(self, values: np.ndarray):
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
                    "tri_s": values[14],
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
            self.basin_structure._Basin.replace(".tif", ".shp"))
        centroid = watershed.centroid
        # Get the epsg code from the shp
        epsg_code = centroid.crs.srs
        x = centroid.x.values.tolist()
        y = centroid.y.values.tolist()
        # Convert the centroid from utm to lat lon
        y, x = projections.utm_to_latlon(y, x, epsg_code)
        # Convert the y to a string to a string and round irt to two decimals
        y_str = str(round(y[0], 2))
        y_str = y_str.replace(".", "")
        # Convert it back to an integer number
        y_int = int(y_str)
        return y_int

    def set_solinitial(self, values: np.ndarray):
        self.solInitial = {"hsini": values[0],
                           "hnini": values[1],
                           "hmini": values[2],
                           "q0": values[3],
                           "tmur": values[4],
                           "tstock": values[5],
                           }

    def set_transfert(self, values: np.ndarray):
        # TODO: Need to find the way to automatically compute the time of concetration
        self._compute_zn()
        self.transfert = {"exxkt": values[0],
                          "zn": values[1]}

    # TODO: This will use the input data from the basin structure.
    def _compute_zn(self, area_th=0.01):
        """_summary_
        This method will compute the concentration time of the given basin.
        """
        # Open the
        CP_fishnet = gpd.read_file(self.basin_structure.cp_fishnet_name)
        cp_struct_name = os.path.join(self.basin_structure.project_path,
                                      "results",
                                      "carreauxPartiels.csv")
        carreaux_partiels = pd.read_csv(cp_struct_name, index_col=0)
        CP_fishnet["x_c"] = CP_fishnet.centroid.x
        CP_fishnet["y_c"] = CP_fishnet.centroid.y
        # Define the pair of coordinates
        # p1 = pd.DataFrame(columns=["X", "Y"])
        # p2 = pd.DataFrame(columns=["X", "Y"])

        # p1 = np.zeros((len(CP_fishnet), 2), dtype=Point)
        # p1 = p1.copy()

        lines = np.zeros(len(CP_fishnet), dtype=LineString)
        max_area = CP_fishnet["cumulArea"].max()
        idx_small_area = CP_fishnet["cumulArea"] > area_th*max_area
        CP_fishnet.index = CP_fishnet["newCPid"].values
        carreaux_partiels.index = CP_fishnet["newCPid"].values
        df_line = pd.DataFrame(columns=["names", "cumulArea"],
                               data=np.empty(shape=[len(CP_fishnet), 2]),
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
                df_line.loc[i, "names"] = f"CP{i} to CP{cp_aval}"
            else:
                p1 = Point(CP_fishnet.loc[i, "x_c"],
                           CP_fishnet.loc[i, "y_c"])
                p2 = Point(CP_fishnet.loc[cp_aval, "x_c"],
                           CP_fishnet.loc[cp_aval, "y_c"])
                lines[i-1] = LineString([p2, p1])
                df_line.loc[i, "names"] = f"CP{i} to CP{cp_aval}"

        gpdf_line = gpd.GeoDataFrame(df_line,
                                     geometry=lines,
                                     crs=CP_fishnet.geometry.crs)

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
        for i in range(len(outlet_routes_complete)):
            routes = outlet_routes_complete.iloc[i, :]
            sub_gpdf = gpdf_line.loc[routes, :]
            lengths.append(sub_gpdf.geometry.length.sum())
        # Find the maximun length value
        idx_longest_path = lengths.index(max(lengths))
        main_route = outlet_routes_complete.iloc[idx_longest_path, :]
        gpdf_line["main_path"] = 0
        gpdf_line.loc[main_route, "main_path"] = 1
        # Open the outlet routes
        streams_file = os.path.join(self.basin_structure.project_path,
                                    "geographic",
                                    "streams_cequeau.shp")
        gpdf_line.to_file(streams_file)
        # ax.quiver(p1[idx_small_area, 0], p1[idx_small_area, 1],
        #           u_c[idx_small_area], v_c[idx_small_area],
        #           scale_units='xy',
        #           angles='xy',
        #           scale=1)
        return False

    def set_fonte(self, values: np.ndarray, model: int):
        # Parameters for the cequeau model
        # DEGREE-DAY
        if model == 1:
            self.fonte = {"cequeau": {
                "strne_s": values[0],
                "tfc_s": values[1],
                "tfd_s": values[2],
                "tsc_s": values[3],
                "tsd_s": values[4],
                "ttd": values[5],
                "tts_s": values[6],
                "jonei": self.option["jonei"],
                "tmur": self.solInitial["tmur"],
                "tstock": self.solInitial["tstock"]}
            }
        # CEMANAIEGE
        elif model == 2:
            # Here we get the mean altitude from the
            # basin structure
            self.fonte = {"cemaNeige": {
                "Kf": values[0],
                "Tf": values[0],
                "CTg": values[0],
                "theta": values[0],
                "QNBV": values[0],
                "Zmed": np.array(self.basin_structure.carreauxPartiels["altitudeMoy"], dtype=np.float32).tolist()
            }
            }

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
                        "corsol": values[8],
                        # This is true if the simulations starts in the winter
                        "panap": 0,
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
