function fixed_struct = fix_struct(carreux_struct)
    % get the table names
    table_names = fieldnames(carreux_struct);
    [~, n] = size(carreux_struct.(table_names{1}));
    % Loop into the struct names
    for i = 1:1:size(table_names, 1)

        % Check the value of this variables because the lists need to be appended differently
        if strcmp(table_names{i}, "idCPsAmont")
            % Loop into each entry
            for j = 1:1:n
                % get the array's column as list
                list_array = carreux_struct.(table_names{i})(j, :);
                fixed_struct(j).(table_names{i}) = list_array;
            end

        else
            % Loop into each entry
            for j = 1:1:n
                fixed_struct(j).(table_names{i}) = carreux_struct.(table_names{i})(j);
            end

        end

    end

end
