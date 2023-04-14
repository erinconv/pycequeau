function grid = upload_meteo(meteo)
    % TODO: This upload process need to catch exceptions,
    % like for example taking into account that all the variables are correctly provided
    % This is the list of the variables fro the cequeau grid format
    var_names = ["pTot", "tMax", "tMin", "pression", ...
                     "rayonnement", "vitesseVent", "nebulosite"];
    ncID = netcdf.open(meteo);
    % Create the meteopointGrille structure for CEQUEAU
    grid = create_grid();

    for idx = 1:1:size(var_names, 2)
        varID = netcdf.inqVarID(ncID, var_names(idx));
        varData = netcdf.getVar(ncID, varID);
        grid.(var_names(idx)) = varData';
    end

    varID = netcdf.inqVarID(ncID, "pasTemp");
    grid.t = netcdf.getVar(ncID, varID);
    netcdf.close(ncID)
end
