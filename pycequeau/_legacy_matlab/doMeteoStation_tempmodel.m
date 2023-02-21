function [stations, meteoStation]=doMeteoStation_tempmodel(CEgrid,R,metStations)
%DOMETEOSTATION Assembles met station data for watershed.  
%   Uses CEgrid and 'metStations' structure from meteoStationData.m 
%   function (see /meteo folder) to assemble stations and meteoStation
%   sub-structures for CEQUEAU
%
%   Input:  'CEgrid'        - CEgrid raster. Same dimensions as FAC/CAT/DEM rasters.
%           'R'             - Raster worldfile giving coordinate system of CAT/FAT/DEM rasters.  Loaded using worldfileread.m function of Matlab Mapping Toolbox
%           'metStations'   - metStations structure generated by meteoStationData.m or appendMeteoStationData.m or mergeMeteo.m (see 'meteo' folder)
%          
%   Output: 'stations'      - CEQUEAU 'stations' sub-structure containing physical attributes for stations
%           'meteoStation'  - CEQUEAU 'meteoStation' sub-structure containing met station records
%           
%   By Stephen Dugdale, 2015-11-12

%sum CE grid to create two vectors that show where CE changes
CEgridX=sum(CEgrid,1);
CEgridY=flipud(sum(CEgrid,2));

%with this bit of the code, we want to renumber the vectors so that they form the indexes for the (i,j) numbering of the CEs (from 10 upwards)

%first, renumber the x-vector
currentValue=CEgridX(1); %set intial 'currentValue'
j=10; %set increment to 10
xSize=0;
for n=1:numel(CEgridX); %scan along vector one cell at a time 
    if CEgridX(n)==0
       CEgridX(n)=NaN;
       continue
    end
    
    if CEgridX(n)==currentValue; %if the value of the cell is the same as 'currentValue', we must still be in the same (Ith) CE...
       CEgridX(n)=j; %...therefore, set the Nth cell to 'j'...
       %currentValue=CEgridX(n); %...and set the new current value equal to this cell.
    else %if the cell value HAS changed (compared to 'currentValue'), we must have changed cell...
       %if n==1 | n==numel(CEgridX)
       %CEgridX(n)=j;
       %else
       currentValue=CEgridX(n); %...therefore, set a new 'currentValue', based on the value of the new cell
       j=j+1; %up the increment by one, because we must be in a new cell 
       CEgridX(n)=j; %set the value of the new cell to this new increment
       xSize=[xSize;n]; %get x dimension size of CE
       %end
    end
end
CEgridX=CEgridX-(min(CEgridX,[], 'omitnan')-10);

xSize=diff(xSize);
xSize(xSize<(max(xSize).*0.99))=max(xSize);
    
%repeat the exact same process for the y-vector
currentValue=CEgridY(1);
j=10;
ySize=0;
for n=1:numel(CEgridY);
    if CEgridY(n)==0
       CEgridY(n)=NaN;
       continue
    end
    
    if CEgridY(n)==currentValue;
       CEgridY(n)=j;
    else
       %if n==1 | n==numel(CEgridY)
       %CEgridY(n)=j;
       %else
       currentValue=CEgridY(n); 
       j=j+1;
       CEgridY(n)=j;
       ySize=[ySize;n]; %get y dimension size of CE
       %end
    end
end
CEgridY=CEgridY-(min(CEgridY, [], 'omitnan')-10);

ySize=diff(ySize);
ySize(ySize<(max(ySize).*0.99))=max(ySize);

CEgridX=CEgridX';%convert CEgridX to column vector

%lengthen grid vectors to include 1:10 and n:n+10 (normally these are
%outside the limits of the CE grid
if numel(xSize)<10
xSize(numel(xSize)+1:10)=mode(xSize);
if size(xSize,2)>size(xSize,1),xSize=xSize';end
xTemp=[0;cumsum(xSize(1:10))]; %find the point at which the i coordinate increments 
else    
xTemp=[0;cumsum(xSize(1:10))]; %find the point at which the i coordinate increments
end
for n=1:numel(xTemp)-1;
CEgridXNeg(xTemp(n)+1:xTemp(n+1),1)=n; %create the negative (ie. i = 1:10) values
end
CEgridXPos=CEgridXNeg+max(CEgridX)-1; %create the positive values

idx=find(CEgridX>10 & CEgridX<max(CEgridX)); %get index of CEgrid positions
CEgridX(CEgridX==10)=[];
CEgridX(CEgridX==max(CEgridX))=[];
CEgridX(isnan(CEgridX))=[];

if ~isempty(CEgridX)
CEgridX_final=(idx(1)-numel(CEgridXNeg):idx(end)+numel(CEgridXPos))';
else
CEgridX_final=(1-numel(CEgridXNeg):numel(CEgridXPos))';
end
%CEgridX_final=(1-numel(CEgridXNeg):numel(CEgridX)+numel(CEgridXPos))'; %give index values to every point within the grid
CEgridX_final(:,2)=[CEgridXNeg;CEgridX;CEgridXPos]; %merge with the 10:n CE grid)

%do the same for the j vectors
if numel(ySize)<10
ySize(numel(ySize)+1:10)=mode(ySize);
if size(ySize,2)>size(ySize,1),ySize=ySize';end
yTemp=[0;cumsum(ySize(1:10))];
else    
yTemp=[0;cumsum(ySize(1:10))];
end
for n=1:numel(yTemp)-1;
CEgridYNeg(yTemp(n)+1:yTemp(n+1),1)=n;
end
CEgridYPos=CEgridYNeg+max(CEgridY)-1;

%flip all vectors upside down because of the way Matlab treats Y axis in rasters
CEgridY=flipud(CEgridY);
CEgridYNeg=flipud(CEgridYNeg);
CEgridYPos=flipud(CEgridYPos);

idx=find(CEgridY>10 & CEgridY<max(CEgridY)); %get index of CEgrid positions
CEgridY(CEgridY==10)=[];
CEgridY(CEgridY==max(CEgridY))=[];
CEgridY(isnan(CEgridY))=[];

if ~isempty(CEgridY)
CEgridY_final=(idx(1)-numel(CEgridYNeg):idx(end)+numel(CEgridYPos))';
else
CEgridY_final=(1-numel(CEgridYNeg):numel(CEgridYPos))';
end
%CEgridY_final=flipud((1-numel(CEgridYNeg):numel(CEgridY)+numel(CEgridYPos))'); %flip CEgridY upside down because of the way Matlab treats Y axis in rasters
CEgridY_final(:,2)=[CEgridYPos;CEgridY;CEgridYNeg];

%convert met station coordinates to UTM
fn=fieldnames(metStations);

for n=1:size(fn,1);
lat(n,1)=metStations.(fn{n}).LATITUDE;
lon(n,1)=metStations.(fn{n}).LONGITUDE;
end
[X Y]=ll2utm(lat,lon,'19N');
Y=lat;
X=lon;
[row, col]=map2pix(R,X,Y);

%get i,j positions of met. stations
%h = waitbar(0,'Getting (I,J) for CEs...');
for n=1:numel(row) %loop through CEs 
    idx=find(CEgridX_final(:,1)==round(col(n)));
    iMet(n,1)=(CEgridX_final(idx,2)); %find i coordinate of Nth CE
    idx=find(CEgridY_final(:,1)==round(row(n)));
    jMet(n,1)=(CEgridY_final(idx,2)); %find j coordinate of Nth CE
    %waitbar(n / numel(row));%update waitbar
end
%close(h);

%find longest timestep with non-NaN start/finish for met stations
for n=1:size(fn,1);
date=metStations.(fn{n}).DATE;
tMax=metStations.(fn{n}).TMAX; tMax(tMax==-9999)=NaN;
tMin=metStations.(fn{n}).TMIN; tMin(tMin==-9999)=NaN;
prcp=metStations.(fn{n}).PRCP; prcp(prcp==-9999)=NaN;
temp=mean([tMax tMin prcp],2,'omitnan');
idx=find((~isnan(temp)));
startdate(n,1)=date(idx(1));
enddate(n,1)=date(idx(end));
rlength(n,1)=enddate(n,1)-startdate(n,1); %update record lengths for stations
end

startdate=min(startdate); %get start date for non-NaN record
enddate=max(enddate); %get end date for non-NaN record
%T=(startdate:enddate)'; %create t vector for meteoStation structure
T = date;
%loop through CEs and find weather stations contained within
num=1;
pTot=[];
tMin=[];
tMax=[];

%do temperature model data
rayonnement=[];
nebulosite=[];
pression=[];
vitesseVent=[];
neige = [];

for i=1:max(CEgridX_final(:,2));
    for j=1:max(CEgridY_final(:,2));
        idx=find(iMet==i & jMet==j); %find weather stations contained within CE(i,j)
        
        %-----get meteo data if multiple weather stations exist in CE(i,j)
        if numel(idx)>1      
        
        pTot=[pTot,nan([numel(T),1])]; %add new columns to pTot,tMin, tMax, snow, snow depth
        tMin=[tMin,nan([numel(T),1])];
        tMax=[tMax,nan([numel(T),1])];
        rayonnement=[rayonnement,nan([numel(T),1])];
        nebulosite=[nebulosite,nan([numel(T),1])];
        pression=[pression,nan([numel(T),1])];
        vitesseVent=[vitesseVent,nan([numel(T),1])];
        neige = [neige,nan([numel(T),1])];
        
        %assemble physiographic data for met stations
        idTemp=[];
        nomTemp=[];
        iTemp=[];
        jTemp=[];
        altitudeTemp=[];
        for a=1:numel(idx);
        idTemp=[idTemp,strtrim(metStations.(fn{idx(a)}).STATION),' '];
        nomTemp=[nomTemp,strtrim(metStations.(fn{idx(a)}).STATION_NAME),' '];
        iTemp=[iTemp iMet(idx(a))];
        jTemp=[jTemp jMet(idx(a))];
        altitudeTemp=[altitudeTemp metStations.(fn{idx(a)}).ELEVATION];
        end
            
        stations(num).id=strtrim(idTemp); %enter appropriate data into the 'stations' structure
        stations(num).nom=strtrim(nomTemp);
        stations(num).i=mean(iTemp);
        stations(num).j=mean(jTemp);
        stations(num).tp=NaN;
        stations(num).altitude=mean(altitudeTemp);
        
        %merge data from met stations
        %get date vector
        dateTemp=[]; 
        for a=1:numel(idx);
        dateTemp=[dateTemp;metStations.(fn{idx(a)}).DATE];
        end
        dateTemp=(min(dateTemp):max(dateTemp))';
        
        %get weather data based on date vector and merge
        pTotTemp=NaN(numel(dateTemp),numel(idx));
        tMinTemp=NaN(numel(dateTemp),numel(idx));
        tMaxTemp=NaN(numel(dateTemp),numel(idx));
        rayonnementTemp=NaN(numel(dateTemp),numel(idx));
        nebulositeTemp=NaN(numel(dateTemp),numel(idx));
        pressionTemp=NaN(numel(dateTemp),numel(idx));
        vitesseVentTemp=NaN(numel(dateTemp),numel(idx));
        neigeTemp=NaN(numel(dateTemp),numel(idx));
        
        for a=1:numel(idx);
        datidx=find(ismember(dateTemp,metStations.(fn{idx(a)}).DATE));
        pTotTemp(datidx,a)=metStations.(fn{idx(a)}).PRCP;
        tMinTemp(datidx,a)=metStations.(fn{idx(a)}).TMIN;
        tMaxTemp(datidx,a)=metStations.(fn{idx(a)}).TMAX;
        if isfield(metStations.(fn{idx(a)}),'SRAD');rayonnementTemp(datidx,a)=metStations.(fn{idx(a)}).SRAD;end
        if isfield(metStations.(fn{idx(a)}),'NEBU');nebulositeTemp(datidx,a)=metStations.(fn{idx(a)}).NEBU;end
        if isfield(metStations.(fn{idx(a)}),'PRES');pressionTemp(datidx,a)=metStations.(fn{idx(a)}).PRES;end
        if isfield(metStations.(fn{idx(a)}),'SRAD');vitesseVentTemp(datidx,a)=metStations.(fn{idx(a)}).WIND;end
        if isfield(metStations.(fn{idx(a)}),'SNOW');neigeTemp(datidx,a)=metStations.(fn{idx(a)}).SNOW;end
        end      
        
        pTotTemp(pTotTemp==-9999)=NaN;
        tMinTemp(tMinTemp==-9999)=NaN;
        tMaxTemp(tMaxTemp==-9999)=NaN;
        rayonnementTemp(rayonnementTemp==-9999)=NaN;
        nebulositeTemp(nebulositeTemp==-9999)=NaN;
        pressionTemp(pressionTemp==-9999)=NaN;
        vitesseVentTemp(vitesseVentTemp==-9999)=NaN; 
        neigeTemp(neigeTemp==-9999)=NaN; 
        
        pTotTemp=mean(pTotTemp,2,'omitnan');
        tMinTemp=mean(tMinTemp,2,'omitnan');
        tMaxTemp=mean(tMaxTemp,2,'omitnan');
        rayonnementTemp=mean(rayonnementTemp,2,'omitnan');
        nebulositeTemp=mean(nebulositeTemp,2,'omitnan');
        pressionTemp=mean(pressionTemp,2,'omitnan');
        vitesseVentTemp=mean(vitesseVentTemp,2,'omitnan');
        neigeTemp=mean(neigeTemp,2,'omitnan');
        
        
        %find where new met data vectors fit within global record
        datidx=find(ismember(T,dateTemp));
        recidx=find(ismember(dateTemp,T));
        
        %add vectors to global record
        pTot(datidx,num)=pTotTemp(recidx); 
        tMin(datidx,num)=tMinTemp(recidx);
        tMax(datidx,num)=tMaxTemp(recidx);
        rayonnement(datidx,num)=rayonnementTemp(recidx); 
        nebulosite(datidx,num)=nebulositeTemp(recidx);
        pression(datidx,num)=pressionTemp(recidx);
        vitesseVent(datidx,num)=vitesseVentTemp(recidx);
        neige(datidx,num)=neigeTemp(recidx);
        
        num=num+1;    
            
        end
        
        
        %-----get weather data if only one station exists within CE(i,j)
        
        if numel(idx)==1
            
        pTot=[pTot,nan([numel(T),1])]; %add new columns to pTot,tMin, tMax, snow, snow depth
        tMin=[tMin,nan([numel(T),1])];
        tMax=[tMax,nan([numel(T),1])];
        rayonnement=[rayonnement,nan([numel(T),1])];
        nebulosite=[nebulosite,nan([numel(T),1])];
        pression=[pression,nan([numel(T),1])];
        vitesseVent=[vitesseVent,nan([numel(T),1])];
        neige=[neige,nan([numel(T),1])];
        
        
        %assemble physiographic data for met stations
        stations(num).id=strtrim(metStations.(fn{idx}).STATION); %enter appropriate data into the 'stations' structure
        stations(num).nom=strtrim(metStations.(fn{idx}).STATION_NAME);
        stations(num).i=iMet(idx);
        stations(num).j=jMet(idx);
        stations(num).tp=NaN;
        stations(num).altitude=metStations.(fn{idx}).ELEVATION;
        datidx=find(ismember(T,metStations.(fn{idx}).DATE));
        recidx=find(ismember(metStations.(fn{idx}).DATE,T));
        
        %create pTot,tMin and tMax vectors
        pTot(datidx,num)=metStations.(fn{idx}).PRCP(recidx); 
        tMin(datidx,num)=metStations.(fn{idx}).TMIN(recidx);
        tMax(datidx,num)=metStations.(fn{idx}).TMAX(recidx);
        
        %do snow vector
        if isfield(metStations.(fn{idx}),'SRAD');rayonnement(datidx,num)=metStations.(fn{idx}).SRAD(recidx);end
        if isfield(metStations.(fn{idx}),'NEBU');nebulosite(datidx,num)=metStations.(fn{idx}).NEBU(recidx);end
        if isfield(metStations.(fn{idx}),'PRES');pression(datidx,num)=metStations.(fn{idx}).PRES(recidx);end
        if isfield(metStations.(fn{idx}),'WIND');vitesseVent(datidx,num)=metStations.(fn{idx}).WIND(recidx);end
        if isfield(metStations.(fn{idx}),'SNOW');neige(datidx,num)=metStations.(fn{idx}).SNOW(recidx);end
        
        num=num+1;
        end
    end   
end

%append data to meteoStation
pTot(pTot==-9999)=NaN; %pTot=pTot./10;
tMin(tMin==-9999)=NaN; %tMin=tMin./10;
tMax(tMax==-9999)=NaN; %tMax=tMax./10;
rayonnement(rayonnement==-9999)=NaN; 
nebulosite(nebulosite==-9999)=NaN;
pression(pression==-9999)=NaN; 
vitesseVent(vitesseVent==-9999)=NaN;
neige(neige==-9999)=NaN;
meteoStation.t=T;
meteoStation.pTot=pTot;
meteoStation.tMin=tMin;
meteoStation.tMax=tMax;
meteoStation.rayonnement=rayonnement;
meteoStation.nebulosite=nebulosite;
meteoStation.pression=pression;
meteoStation.vitesseVent=vitesseVent;
meteoStation.neige=neige;

%do annual precip
tVec=datevec(T); %seperate serial dates into vector
years=unique(tVec(:,1)); %extract years from vector

for n=1:numel(stations);
    annualPrecip=[]; %seed annualPrecip variable
    for m=1:numel(years); %loop through years
        tempPrecip=pTot(tVec(:,1)==years(m),n); %get vector of precip in Mth year
        tempPrecip(isnan(tempPrecip))=[]; %remove NaNs
        if numel(tempPrecip)>=355 %if 355 or more observations exist...
        annualPrecip=[annualPrecip;sum(tempPrecip)]; %...sum the annual precipitation and add to annualPrecip vector   
        end
    end
    stations(n).tp=mean(annualPrecip); %get the mean of annualPrecip vector for Nth station and append to stations structure
end      


end