function meteoPointGrille = doMeteoInterpolation(inputStruct,corr)
    %DOMETEOINTERPOLATION Interpolates met station data
    %   When normal CEQUEAU interpolator won't work because of stations with 
    %   missing records, this function will do a simple three nearest neighbour
    %   interpolation of met station data that will output a meteoPointGrille
    %   file compatible with CEQUEAU.  Input data is CEQUEAU structure 
    %   containing stations and meteoStation sub-structures.  Output is 
    %   meteoPointGrille of met. data.
    %
    %   meteoPointGrille = doMeteoInterpolation(inputStruct,corr)
    %
    %   Input:  'inputStruct'       -   CEQUEAU input structure containing meteo data).
    %           'corr'              -   switch whether we decide to use altitude correction for temperature data (based on coet/coep values in inputStruct.parametres.interpolation).   if 'corr' == 1, then correction will be applied.
    %
    %   Output: 'meteoPointGrille'  -   meteoPointGrille file compatible with CEQUEAU  
    %
    %   By Steve Dugdale, 2015.  Based on function 'interpolerNearest_3' by Marco Latraverse
    
    if nargin<2
       corr=0;
    end
    
    %get CE data
    altCE=[inputStruct.bassinVersant.carreauxEntiers.altitude];
    iCE=[inputStruct.bassinVersant.carreauxEntiers.i]';
    jCE=[inputStruct.bassinVersant.carreauxEntiers.j]';
    
    %get station data
    iStation=[inputStruct.stations.i]';
    jStation=[inputStruct.stations.j]';
    tpStation=[inputStruct.stations.tp];
    altStation=[inputStruct.stations.altitude];
    
    %get meteo records
    pTot=inputStruct.meteoStation.pTot;
    tMin=inputStruct.meteoStation.tMin;
    tMax=inputStruct.meteoStation.tMax;
    rayonnement = inputStruct.meteoStation.rayonnement;
    nebulosite = inputStruct.meteoStation.nebulosite;
    pression =inputStruct.meteoStation.pression;
    vitesseVent = inputStruct.meteoStation.vitesseVent;
    %neige = inputStruct.meteoStation.neige;
    
    %get correction factor
    coep=inputStruct.parametres.interpolation.coep;
    coet=inputStruct.parametres.interpolation.coet;
    
    for n=1:size(tMax,1);
    
    %---do tMax---
    tMax_timestep=tMax(n,:)'; %get tMax values at Nth timestep
    if size(tMax_timestep(tMax_timestep==NaN),1) > 0
    %if anyEq(tMax_timestep,NaN)
    %'space'
    end
    
    i_timestep_tMax=iStation(~isnan(tMax_timestep)); %get i coordinates of stations with non-NaN values
    j_timestep_tMax=jStation(~isnan(tMax_timestep)); %get j coordinates of stations with non-NaN values
    tMax_timestep(isnan(tMax_timestep))=[]; %remove NaN values
    
    [D,Idistance] = pdist2([i_timestep_tMax j_timestep_tMax], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE and get distance
    
    D = D + .00001; %add tiny amount to get rid of 'divide by zero' errors
    
    W = (1./D) ./ repmat(sum(1./D),3,1); %compute weighting for each station based on distances
    
    if corr==1
    FCt = coet.*(repmat(altCE,3,1)-altStation(Idistance))/1000; %compute temperature correction factor for altitude
    meteoPointGrille.tMax(n,:) = sum(W .* (FCt+tMax_timestep(Idistance))); %apply correction factor and interpolate tMAX for each CE
    else
    meteoPointGrille.tMax(n,:) = sum(W .* tMax_timestep(Idistance)); %apply correction factor and interpolate tMAX for each CE   
    end
    
    %---do tMin---
    
    tMin_timestep=tMin(n,:)';
    i_timestep_tMin=iStation(~isnan(tMin_timestep));
    j_timestep_tMin=jStation(~isnan(tMin_timestep));
    tMin_timestep(isnan(tMin_timestep))=[];
    
    [D,Idistance] = pdist2([i_timestep_tMin j_timestep_tMin], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    D = D + .00001;
    W = (1./D) ./ repmat(sum(1./D),3,1);
    
    if corr==1
    FCt = coet.*(repmat(altCE,3,1)-altStation(Idistance))/1000;
    meteoPointGrille.tMin(n,:) = sum(W .* (FCt+tMin_timestep(Idistance)));
    else
    meteoPointGrille.tMin(n,:) = sum(W .* tMin_timestep(Idistance));   
    end
    
    %---do pTot---
    
    pTot_timestep=pTot(n,:)';
    i_timestep_pTot=iStation(~isnan(pTot_timestep));
    j_timestep_pTot=jStation(~isnan(pTot_timestep));
    pTot_timestep(isnan(pTot_timestep))=[];
    
    [D,Idistance] = pdist2([i_timestep_pTot j_timestep_pTot], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    D = D + .00001;
    W = (1./D) ./ repmat(sum(1./D),3,1);
    
    meteoPointGrille.pTot(n,:) = sum(W .* pTot_timestep(Idistance));   
    
    %Solar radiation
    rayonnement_timestep=rayonnement(n,:)';
    i_timestep_rayonnement=iStation(~isnan(rayonnement_timestep));
    j_timestep_rayonnement=jStation(~isnan(rayonnement_timestep));
    rayonnement_timestep(isnan(rayonnement_timestep))=[];
    
    [D,Idistance] = pdist2([i_timestep_rayonnement j_timestep_rayonnement], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    D = D + .00001;
    W = (1./D) ./ repmat(sum(1./D),3,1);
    
    meteoPointGrille.rayonnement(n,:) = sum(W .* rayonnement_timestep(Idistance));   
    
    %Cloud covering
    nebulosite_timestep=nebulosite(n,:)';
    i_timestep_nebulosite=iStation(~isnan(nebulosite_timestep));
    j_timestep_nebulosite=jStation(~isnan(nebulosite_timestep));
    nebulosite_timestep(isnan(nebulosite_timestep))=[];
    
    [D,Idistance] = pdist2([i_timestep_nebulosite j_timestep_nebulosite], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    D = D + .00001;
    W = (1./D) ./ repmat(sum(1./D),3,1);
    
    meteoPointGrille.nebulosite(n,:) = sum(W .* nebulosite_timestep(Idistance));   
    
    
    %Pressure
    pression_timestep=pression(n,:)';
    i_timestep_pression=iStation(~isnan(pression_timestep));
    j_timestep_pression=jStation(~isnan(pression_timestep));
    pression_timestep(isnan(pression_timestep))=[];
    
    [D,Idistance] = pdist2([i_timestep_pression j_timestep_pression], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    D = D + .00001;
    W = (1./D) ./ repmat(sum(1./D),3,1);
    
    meteoPointGrille.pression(n,:) = sum(W .* pression_timestep(Idistance));   
    
    
    %Wind
    vitesseVent_timestep=vitesseVent(n,:)';
    i_timestep_vitesseVent=iStation(~isnan(vitesseVent_timestep));
    j_timestep_vitesseVent=jStation(~isnan(vitesseVent_timestep));
    vitesseVent_timestep(isnan(vitesseVent_timestep))=[];
    
    [D,Idistance] = pdist2([i_timestep_vitesseVent j_timestep_vitesseVent], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    D = D + .00001;
    W = (1./D) ./ repmat(sum(1./D),3,1);
    
    meteoPointGrille.vitesseVent(n,:) = sum(W .* vitesseVent_timestep(Idistance));
    
    
    %Snow
    %neige_timestep=neige(n,:)';
    %i_timestep_neige=iStation(~isnan(neige_timestep));
    %j_timestep_neige=jStation(~isnan(neige_timestep));
    %neige_timestep(isnan(neige_timestep))=[];
    
    %[D,Idistance] = pdist2([i_timestep_neige j_timestep_neige], [iCE jCE], 'euclidean', 'Smallest', 3); %find nearest three stations to each CE
    %D = D + .00001;
    %W = (1./D) ./ repmat(sum(1./D),3,1);
    
    %meteoPointGrille.neige(n,:) = sum(W .* neige_timestep(Idistance));
    
    end