import numpy as np


def send_values_test() -> tuple:
    """_summary_

    Returns:
        tuple: _description_
    """

    flow_parameters = [
        # CIN Infiltration coefficient from SOL (upper) reservoir to the NAPPE (lower) reservoir
        0.4,
        # CVMAR Drainage coefficient for the LACS  & MARAIS (lakes and marshes) reservoir
        0.5036,
        # CVNB Lower Drainage coefficient for the NAPPE (lower) reservoir
        0.0156,
        # CVNH Upper drainage coefficient for the NAPPE (lower) reservoir
        0.17382,
        # CVSB Lower Drainage coefficient for the SOL (upper) reservoir
        0.102233,
        # CVSI Upper drainage coefficient for the SOL (upper) reservoir
        0.44716,
        # XINFMA Maximum allowable daily infiltration from SOL to NAPPE (mm/day)
        9.5775,
        # HINF Infiltration threshold (minimum water level) from SOL to NAPPE (mm)
        11.718,
        # HINT Intermediate level for SOL reservoir drainage (mm)
        107.4,
        # HMAR Drainage level threshold for LACS et MARAIS (lakes and marshes) reservoir (mm)
        370.77,
        # HNAP Upper drainage level threshold for the NAPPE reservoir (mm)
        94.943,
        # HPOT Threshold of minimum water level to allow water evapotranspiration a the potential rate  (mm)
        77.008,
        # HSOL Height of reservoir SOL (mm)
        181.18,
        # HRIMP Minimum water level required to initiate runoff on impervious surfaces (mm)
        5.0666,
        # TRI percentage of impermeable surface
        0
    ]

    snow_parameters = [
        # STRNE Snow-rain temperature threshold (°C)
        -0.8285,
        # TFC Potential melting rate in forest  (mm/°C/jour)
        8.6557,
        # TFD Potential melting rate in open (no canopy) areas (mm/°C/jour)
        6.3987,
        # TSC Minimum temperature threshold to initiate snowmelt in forest (°C)
        -1.8379,
        # TSD Minimum temperatue threshold to initiate snowmelt in open areas (°C)
        -0.22692,
        # TTD Heat deficit coefficient (°C)
        2.8166,
        # TTS Minimum temperature for snow stock ripening (°C)
        0.45228
    ]
    evapo_parameters = [
        # EVNAP Fraction of evapotranspiration taken for the NAPPE reservoir
        0.1481,
        # XAA Thorntwaite exponent
        0.9017,
        # XIT Thorntwaite Index
        9.861,
    ]

    initial_conditions = [
        # HSNI
        5.0,
        # HNINI
        5.0,
        # HMINI
        100.0,
        # q0
        10.0,
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
        0.029266,
        # ZN Time of concentration of the basin. Can be calculated also
        4.641066925
    ]

    temperature_params = [
        # COPROM  Coefficient defining minimum river depth as a ratio of width
        2,
        # COLARG Coefficient defining minimum river width
        1.39204613161842,
        # CRAYSO Weighting coefficient for solar (short wave) radiation
        2.22370022406346,
        # CRAYIN Weighting coefficient for infrared radiation
        1.52279743689713,
        # CEVAPO Weighting coefficient for evaporation (latent heat)
        0.5,
        # CCONVE Weighting coefficient for convection (sensible heat)
        1.63455345799646,
        # CRIGEL Freeze criterion for all whole squares (minimum amount of snow in mm)
        74.3622181298935,
        # TNAP Groundwater temperature (°C)
        7.80799635589295,
        # BASSOL Total precipitation required to dectect days with low solar radiation (mm)
        9.13537208064851,
        # Correction du rayonnement solaire moyen (RSM) pour les jours sanspluie (RSM&(1+CORSOL)) et les jours de fortes pluies(RSM&(1-CORSOL)) (varie entre 0,0 et 1,0).
        0.0494176705022366
    ]
    # Convert the parameter list into the required format for the CEQUEAU model
    flow_parameters = np.array(flow_parameters, dtype=np.double).tolist()
    evapo_parameters = np.array(evapo_parameters, dtype=np.double).tolist()
    initial_conditions = np.array(
        initial_conditions, dtype=np.double).tolist()
    snow_parameters = np.array(snow_parameters, dtype=np.double).tolist()
    transferts = np.array(transferts, dtype=np.double).tolist()
    temperature_params = np.array(
        temperature_params, dtype=np.double).tolist()
    simulation_options = np.array(simulation_options, dtype=np.int8).tolist()
    return flow_parameters, evapo_parameters, initial_conditions, snow_parameters, simulation_options, transferts, temperature_params
