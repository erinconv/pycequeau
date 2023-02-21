function visualiseMeteo(inputStruct, type)
    %VISUALISEMETEO Provides animation of met data for each timestep
    %   Draws and animates met data from CEQUEAU meteoPointGrille structure.
    %   One day of weather = 1/100th of a second of animation.  Each pixel of
    %   animation is one CE in watershed
    %
    %   visualiseMeteo(inputStruct,type)
    %
    %   Input:  'inputStruct'   - CEQUEAU structure.  MUST contain meteoPointGrille sub-structure
    %           'type'          - Weather type.  Either tMin, tMax or pTot
    %
    %   By Stephen Dugdale, 2015-04-24

    if (strcmp(type, 'tMin') | strcmp(type, 'tMax') | strcmp(type, 'pTot') | strcmp(type, 'vitesseVent') | strcmp(type, 'pression') | strcmp(type, 'nebulosite') | strcmp(type, 'rayonnement')) == 0
        warndlg('Unrecognised data type', 'Warning');
        return
    end

    i = [inputStruct.bassinVersant.carreauxEntiers.i]; %get i coords of CEs
    j = [inputStruct.bassinVersant.carreauxEntiers.j]; %get j coords of CEs

    if min(i) > 1
        i = i - (min(i) - 1);
    end

    if min(j) > 1
        j = j - (min(j) - 1);
    end

    startdate = find(inputStruct.meteoStation.t == inputStruct.execution.dateDebut);
    enddate = find(inputStruct.meteoStation.t == inputStruct.execution.dateFin);

    minVal = min(min(inputStruct.meteoPointGrille.(type)));
    maxVal = max(max(inputStruct.meteoPointGrille.(type)));

    figure

    for m = startdate:enddate; %change the '2940' back to 1 ...
            metgrid = nan(max(j), max(i));

        for n = 1:size(inputStruct.meteoPointGrille.(type), 2);
            metgrid(j(n), i(n)) = inputStruct.meteoPointGrille.(type)(m, n);
        end

        h1 = imagesc(metgrid); colormap(jet);
        axis equal
        %if strcmp(type,'tMin') | strcmp(type,'tMax')
        %caxis([-35 35]);
        %else
        %caxis([0 200]);
        %end
        caxis([minVal maxVal]);
        h2 = colorbar;
        set(gca, 'ydir', 'normal');
        h3 = text(2, 2, datestr(inputStruct.meteoStation.t(m), 'yyyy-mm-dd'), 'color', 'w');
        drawnow;
        pause(0.001);
    end
