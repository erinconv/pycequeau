function export_results_CE(resu, path, decade)
    % Create the netcdf
    name = strcat("CE_results", decade, ".nc");
    file_name = fullfile(path, name);
    ncID = netcdf.create(file_name, "NETCDF4");
    % Define the dimensions
    Temp_dim = netcdf.defDim(ncID, "pasTemp", size(resu.t, 1));
    Carreux_dim = netcdf.defDim(ncID, "CEid", size(resu.CE.CEs, 2));
    varid = netcdf.defVar(ncID, "pasTemp", 'NC_FLOAT', Temp_dim);
    netcdf.putVar(ncID, varid, resu.t)
    varid = netcdf.defVar(ncID, "CEid", 'NC_INT', Carreux_dim);
    netcdf.putVar(ncID, varid, int16(resu.CE.CEs))
    CE_names = ["Sol", "Nappe", "LacsMarais", "stockNeigeForet", "evapo"];

    for idx = 1:1:size(CE_names, 2)
        varid = netcdf.defVar(ncID, CE_names(idx), 'NC_FLOAT', [Carreux_dim, Temp_dim]);
        netcdf.putVar(ncID, varid, resu.CE.(CE_names(idx))')
    end

    netcdf.close(ncID)
end
