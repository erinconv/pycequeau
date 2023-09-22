import numpy as np


def send_values_test() -> tuple:
    """_summary_

    Returns:
        tuple: _description_
    """

    flow_parameters = [
        # CIN Infiltration coefficient from SOL (upper) reservoir to the NAPPE (lower) reservoir
        0.6846,
        # CVMAR Drainage coefficient for the LACS  & MARAIS (lakes and marshes) reservoir
        0.997,
        # CVNB Lower Drainage coefficient for the NAPPE (lower) reservoir
        0.1128,
        # CVNH Upper drainage coefficient for the NAPPE (lower) reservoir
        0.011,
        # CVSB Lower Drainage coefficient for the SOL (upper) reservoir
        0.01,
        # CVSI Upper drainage coefficient for the SOL (upper) reservoir
        0.2256,
        # XINFMA Maximum allowable daily infiltration from SOL to NAPPE (mm/day)
        40,
        # HINF Infiltration threshold (minimum water level) from SOL to NAPPE (mm)
        85.491,
        # HINT Intermediate level for SOL reservoir drainage (mm)
        40.035,
        # HMAR Drainage level threshold for LACS et MARAIS (lakes and marshes) reservoir (mm)
        346.91,
        # HNAP Upper drainage level threshold for the NAPPE reservoir (mm)
        100.02,
        # HPOT Threshold of minimum water level to allow water evapotranspiration a the potential rate  (mm)
        110.62,
        # HSOL Height of reservoir SOL (mm)
        56.75,
        # HRIMP Minimum water level required to initiate runoff on impervious surfaces (mm)
        10,
        # TRI percentage of impermeable surface
        0
    ]

    snow_parameters = [
        # STRNE Snow-rain temperature threshold (°C)
        -0.17,
        # TFC Potential melting rate in forest  (mm/°C/jour)
        4,
        # TFD Potential melting rate in open (no canopy) areas (mm/°C/jour)
        4.76,
        # TSC Minimum temperature threshold to initiate snowmelt in forest (°C)
        -1.4,
        # TSD Minimum temperatue threshold to initiate snowmelt in open areas (°C)
        -0.18,
        # TTD Heat deficit coefficient (°C)
        -1.06,
        # TTS Minimum temperature for snow stock ripening (°C)
        -2.54
    ]
    evapo_parameters = [
        # EVNAP Fraction of evapotranspiration taken for the NAPPE reservoir
        0.3925,
        # XAA Thorntwaite exponent
        0.846,
        # XIT Thorntwaite Index
        8.1453,
    ]

    initial_conditions = [
        # HSNI
        200,
        # HNINI
        110,
        # HMINI
        350,
        # q0
        10,
        # TMUR
        0,
        # TSTOCK
        0,
    ]

    simulation_options = [
        # MODULEFONTE
        1,
        # MODULEEVAPO
        1,
        # CALCULQUALITE
        0
    ]

    transferts = [
        # EXXKT Transfer coefficient from one partial square to another
        0.001,
        # ZN Time of concentration of the basin. Can be calculated also
        4.641066925
    ]

    temperature_params = [
        # COPROM  Coefficient defining minimum river depth as a ratio of width
        4,
        # COLARG Coefficient defining minimum river width
        3.5,
        # CRAYSO Weighting coefficient for solar (short wave) radiation
        3.2915,
        # CRAYIN Weighting coefficient for infrared radiation
        0.461,
        # CEVAPO Weighting coefficient for evaporation (latent heat)
        1.1774,
        # CCONVE Weighting coefficient for convection (sensible heat)
        2.5407,
        # CRIGEL Freeze criterion for all whole squares (minimum amount of snow in mm)
        896.5066,
        # TNAP Groundwater temperature (°C)
        8.8115,
        # BASSOL Total precipitation required to dectect days with low solar radiation (mm)
        17.4304,
        # Correction du rayonnement solaire moyen (RSM) pour les jours sanspluie (RSM&(1+CORSOL)) et les jours de fortes pluies(RSM&(1-CORSOL)) (varie entre 0,0 et 1,0).
        0
    ]
    # Convert the parameter list into the required format for the CEQUEAU model
    flow_parameters = np.array(flow_parameters, dtype=np.float32).tolist()
    evapo_parameters = np.array(evapo_parameters, dtype=np.float32).tolist()
    initial_conditions = np.array(
        initial_conditions, dtype=np.float32).tolist()
    snow_parameters = np.array(snow_parameters, dtype=np.float32).tolist()
    transferts = np.array(transferts, dtype=np.float32).tolist()
    temperature_params = np.array(
        temperature_params, dtype=np.float32).tolist()
    simulation_options = np.array(simulation_options, dtype=np.int8).tolist()
    return flow_parameters, evapo_parameters, initial_conditions, snow_parameters, simulation_options, transferts, temperature_params


send_values_test()
