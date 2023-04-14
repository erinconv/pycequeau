% add paths to the functions and dependencies
addpath("functions\")
addpath("jsonlab-2.0\")
% Load the basin file
load('Meleze-1979-2020.mat')
project_path = '\home\erinconv\01-PhD\River';
% Import json parameters data
params_json = fullfile(project_path, 'results', 'parameters.json');
params_struct = loadjson(params_json);
% Import json bassinversant data
bassin_json = fullfile(project_path, 'results', 'bassinVersant.json');
bassin_struct = loadjson(bassin_json);
% fix the table in both, the carreuxpartiels and carreuxentiers structure
bassin_struct.carreauxEntiers = fix_struct(bassin_struct.carreauxEntiers);
bassin_struct.carreauxPartiels = fix_struct(bassin_struct.carreauxPartiels);
% define the start and end dates
execution.dateDebut = datenum(1979, 01, 01);
execution.dateFin = datenum(2020, 12, 31);
% Choose whether you want or not to simulate water temperature
% no = 0; yes = 1
params_struct.option.calculQualite = 0;
% Upload the meteo files
meteo_file = fullfile(project_path,'meteo','meteo_cequeau.nc');
meteo_grid = upload_meteo(meteo_file);
%% Run simulations
% Add path to the folder where the CEQUEAU model binary file is stored
addpath("01-CEQUEAU\")
[y.etatsCE, y.etatsCP, y.etatsFonte, y.etatsEvapo, y.etatsBarrage, y.pasDeTemps, ...
     y.avantassimilationssCE, y.avantassimilationssFonte, ...
     y.avantassimilationssEvapo, y.etatsQualCP, y.avAssimQual] = ...
    cequeauQuantiteMex_v461(execution, params_struct,bassin_struct, meteo_grid, [], []);