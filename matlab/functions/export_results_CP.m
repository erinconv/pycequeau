function export_results_CP(resu, path, decade)
    % Create the netcdf
    name = strcat("CP_results", decade, ".nc");
    file_name = fullfile(path, name);
    ncID = netcdf.create(file_name, "NETCDF4");
    % Define the dimensions
    Temp_dim = netcdf.defDim(ncID, "pasTemp", size(resu.t, 1));
    Carreux_dim = netcdf.defDim(ncID, "CEid", size(resu.CP.CPs, 2));
    varid = netcdf.defVar(ncID, "pasTemp", 'NC_FLOAT', Temp_dim);
    netcdf.putVar(ncID, varid, resu.t)
    varid = netcdf.defVar(ncID, "CEid", 'NC_INT', Carreux_dim);
    netcdf.putVar(ncID, varid, int16(resu.CP.CPs))
    %     resu.CP.ruiss, resu.CP.nappe, resu.CP.hypo, resu.CP.lacma, resu.CP.radso, resu.CP.radin, resu.CP.evap, resu.CP.conv
    %     CP_names = ["debit","temperature","volume","ruiss","nappe","hypo","lacma","radso","radin","evap","conv"];
    CP_names = ["debit", "temperature", "volume", "ruiss", "nappe", "hypo", "lacma", "radso", "radin", "evap", "conv"];

    for idx = 1:1:size(CP_names, 2)
        varid = netcdf.defVar(ncID, CP_names(idx), 'NC_FLOAT', [Carreux_dim, Temp_dim]);
        netcdf.putVar(ncID, varid, resu.CP.(CP_names(idx))')
    end

    netcdf.close(ncID)
end
