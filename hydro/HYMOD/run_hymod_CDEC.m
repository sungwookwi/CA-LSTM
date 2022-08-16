
% HYMOD Run for CALFEWS Watersheds

clear; clc
addpath('G:\CA-LSTM\HYMOD\hymod_module')

trial_all = {'trial_1','trial_2','trial_3','trial_4','trial_5','trial_6','trial_7','trial_8','trial_9','trial_10'};

for ttt = 1%9:10
    
    trial_sel = trial_all{ttt};

% tic;
%% HRU Information File for the whole CALFEWS area
hruinfo_calfews = load('G:\CA-LSTM\HYMOD\hruinfo\HRUinfo_ALL.txt');
hru_lat_calfews = hruinfo_calfews(:,1); % HRU Lat
hru_lon_calfews = hruinfo_calfews(:,2); % HRU Long


%% Calibration result file to load optimal parameters
fid = fopen(['G:\CA-LSTM\HYMOD\calibration\',trial_sel,'\hymod_ga_cdec_pool_optpar.txt']);
hru_par_calfews = cell(length(hru_lat_calfews),17);
while ~feof(fid)
    str = fgets(fid);
    if contains(str,'HYMOD OPTIMAL PARAMETERS')
        fgets(fid);
        fgets(fid);
        n=1;
        for jj = 1:length(hru_lat_calfews)
            str = fgets(fid);
            hru_par_calfews(n,:) = textscan(str,['%s %s '...
                '%f %f %f %f %f %f %f %f %f %f '...
                '%f %f %f %f %f']);
            n=n+1;
        end
        break;
    end
end
fclose(fid);


%% HRU Information File for each watershed
cdec_id = {'SHA','BND','ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB'};
for i=1:length(cdec_id)
    hruinfo = load(['G:\CA-LSTM\HYMOD\hruinfo\HRUinfo_',cdec_id{i},'.txt']);
    % HRU Lat
    eval(['hru_lat_',cdec_id{i},' = hruinfo(:,1);']) 
    % HRU Long
    eval(['hru_lon_',cdec_id{i},' = hruinfo(:,2);']) 
    % HRU Area fraction
    eval(['hru_area_',cdec_id{i},' = hruinfo(:,3);']) 
    % HRU flow length to the basin outlet (m)
    eval(['hru_flen_',cdec_id{i},' = hruinfo(:,4);']) 
end


%% Run SACSMA for the entire CALFEWS watersheds
% Simulation period
sim_startdate = [1987 10 1];
sim_enddate = [2013 9 30];
sim_datemat = datevec(datenum(sim_startdate):datenum(sim_enddate));

% Module initial state
IniState = [0,0,100];

% Sacrameto HRU output initialization
hru_surf = zeros(size(sim_datemat,1),length(hru_lat_calfews));
hru_base = zeros(size(sim_datemat,1),length(hru_lat_calfews));
hru_aet = zeros(size(sim_datemat,1),length(hru_lat_calfews));
hru_swe = zeros(size(sim_datemat,1),length(hru_lat_calfews));
tic;
% parfor i = 1:length(hru_lat_calfews)
for i = 1:length(hru_lat_calfews)    
    hru_meteo = load(['E:\StochasticWatershedModeling\VIC\hru_meteo\hymod\meteo_',num2str(hru_lat_calfews(i),'%1.6f'),'_',num2str(hru_lon_calfews(i),'%1.6f')]); % yr mon day pr(mm) tas(C)
    sind = find(hru_meteo(:,1)==sim_startdate(1)&hru_meteo(:,2)==sim_startdate(2)&hru_meteo(:,3)==sim_startdate(3));
    eind = find(hru_meteo(:,1)==sim_enddate(1)&hru_meteo(:,2)==sim_enddate(2)&hru_meteo(:,3)==sim_enddate(3));
    pr = hru_meteo(sind:eind,4); % HRU pr for the simulation period
    tas = hru_meteo(sind:eind,5); % HRU tas for the simulation period
    
    selind = i;
    % Module parameters
    Coeff = hru_par_calfews{selind,3};
    SMA_Par = [hru_par_calfews{selind,4},hru_par_calfews{selind,5},hru_par_calfews{selind,6},hru_par_calfews{selind,7},hru_par_calfews{selind,8},hru_par_calfews{selind,9},hru_par_calfews{selind,10}];
    SnowPar = [hru_par_calfews{selind,11},hru_par_calfews{selind,12},hru_par_calfews{selind,13}];
    
    % Run modules: Hamon -> Snow -> HYMOD
    pet = pet_hamon_mex(sim_startdate, sim_enddate, tas, hru_lat_calfews(selind), Coeff);  
    [snow_outflow, snow_swe, snow_melt] = snow_ddf1_mex(pr, tas, SnowPar, IniState(2));
    [surf, base, aet] = sma_hymod2(snow_outflow, pet, snow_swe, SMA_Par, [IniState(1),IniState(3)]);
    
    hru_surf(:,i) = surf;
    hru_base(:,i) = base;
    hru_aet(:,i) = aet;
	hru_swe(:,i) = snow_swe;
 
end
toc;

%% Run Lohmann Routing for each watershed
tic;
for k = 1:length(cdec_id)

    eval(['hru_flen_sta = hru_flen_',cdec_id{k},';'])
    eval(['hru_lat_sta = hru_lat_',cdec_id{k},';'])
    eval(['hru_lon_sta = hru_lon_',cdec_id{k},';'])
    eval(['hru_area_sta = hru_area_',cdec_id{k},';'])
    
    hru_flow_sta = zeros(size(sim_datemat,1),length(hru_lat_sta));
    hru_aet_sta = zeros(size(sim_datemat,1),length(hru_lat_sta));
	hru_swe_sta = zeros(size(sim_datemat,1),length(hru_lat_sta));

    [~,minind] = min(hru_flen_sta);
    isOutlet = zeros(length(hru_flen_sta),1);
    isOutlet(minind) = 1;

    %parfor i = 1:length(hru_lat_sta)
    for i = 1:length(hru_lat_sta)
        selind = 0;
        for j=1:length(hru_lat_calfews)
            if strcmp(hru_par_calfews{j,1}{1},num2str(hru_lat_sta(i),'%2.6f')) && strcmp(hru_par_calfews{j,2}{1},num2str(hru_lon_sta(i),'%2.6f'))
                selind = j;
                break;
            end
        end
        route_par = [hru_par_calfews{selind,14},hru_par_calfews{selind,15},hru_par_calfews{selind,16},hru_par_calfews{selind,17}];
        [totflow, baseflow] = rout_lohmann_mex(hru_surf(:,selind), hru_base(:,selind), hru_flen_sta(i), route_par, isOutlet(i));
        hru_flow_sta(:,i) = totflow * hru_area_sta(i)/sum(hru_area_sta);
        hru_aet_sta(:,i) = hru_aet(:,selind) * hru_area_sta(i)/sum(hru_area_sta);
		hru_swe_sta(:,i) = hru_swe(:,selind) * hru_area_sta(i)/sum(hru_area_sta);
    end

    eval(['simflow_',cdec_id{k},' = sum(hru_flow_sta,2);'])
    eval(['simaet_',cdec_id{k},' = sum(hru_aet_sta,2);'])
	eval(['simswe_',cdec_id{k},' = sum(hru_swe_sta,2);'])
    
end
toc; % takes about 680 sec

%% Validation
cal_startdate_other = [1988 10 1];
cal_startdate_BND = [1999 10 1];
cal_enddate = [2003 9 30];

val_startdate = [2003 10 1];
val_enddate = [2013 9 30];

kge_nse_cal = nan(length(cdec_id),2);
kge_nse_val = nan(length(cdec_id),2);
for i=1:length(cdec_id)
    
    eval(['simflow_sta = simflow_',cdec_id{i},';'])  
    eval(['simaet_sta = simaet_',cdec_id{i},';'])
	eval(['simswe_sta = simswe_',cdec_id{i},';'])
    
    if strcmp(cdec_id{i},'BND')
        cal_startdate = cal_startdate_BND;
    else
        cal_startdate = cal_startdate_other;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cal_datemat_cal = datevec(datenum(cal_startdate):datenum(cal_enddate));
    
    sind = find(sim_datemat(:,1)==cal_startdate(1)&sim_datemat(:,2)==cal_startdate(2)&sim_datemat(:,3)==cal_startdate(3));
    eind = find(sim_datemat(:,1)==cal_enddate(1)&sim_datemat(:,2)==cal_enddate(2)&sim_datemat(:,3)==cal_enddate(3));
    simflow_sta_cal = simflow_sta(sind:eind); % Flow simulation for the validation period
    simaet_sta_cal = simaet_sta(sind:eind);
    
    % CDEC FNF for FOL
    cdec_fnf_sta = load(['G:\CA-LSTM\CDEC\FNF_mm_QC\FNF_',cdec_id{i},'_mm.txt']); % [yr mon day flow(mm)]
    sind = find(cdec_fnf_sta(:,1)==cal_startdate(1)&cdec_fnf_sta(:,2)==cal_startdate(2)&cdec_fnf_sta(:,3)==cal_startdate(3));
    eind = find(cdec_fnf_sta(:,1)==cal_enddate(1)&cdec_fnf_sta(:,2)==cal_enddate(2)&cdec_fnf_sta(:,3)==cal_enddate(3));
    cdec_fnf_sta_cal = cdec_fnf_sta(sind:eind,4);
    
    % NSE: any negative FLOW_cal_obs values (either missing or FNF calculation) are excluded
    nse = 1 - (mean((simflow_sta_cal(cdec_fnf_sta_cal>=0)-cdec_fnf_sta_cal(cdec_fnf_sta_cal>=0)).^2)/var(cdec_fnf_sta_cal(cdec_fnf_sta_cal>=0),1));
    
    % KGE
    mean_ratio = mean(simflow_sta_cal(cdec_fnf_sta_cal>=0))/mean(cdec_fnf_sta_cal(cdec_fnf_sta_cal>=0));
    std_ratio = std(simflow_sta_cal(cdec_fnf_sta_cal>=0),1)/std(cdec_fnf_sta_cal(cdec_fnf_sta_cal>=0),1);
    lin_corr = corr(simflow_sta_cal(cdec_fnf_sta_cal>=0),cdec_fnf_sta_cal(cdec_fnf_sta_cal>=0));
    kge = 1-((1-mean_ratio)^2+(1-std_ratio)^2+(1-lin_corr)^2)^0.5;
        
    kge_nse_cal(i,1) = kge;
    kge_nse_cal(i,2) = nse;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cal_datemat_val = datevec(datenum(val_startdate):datenum(val_enddate));
    
    sind = find(sim_datemat(:,1)==val_startdate(1)&sim_datemat(:,2)==val_startdate(2)&sim_datemat(:,3)==val_startdate(3));
    eind = find(sim_datemat(:,1)==val_enddate(1)&sim_datemat(:,2)==val_enddate(2)&sim_datemat(:,3)==val_enddate(3));
    simflow_sta_val = simflow_sta(sind:eind); % Flow simulation for the validation period
    simaet_sta_val = simaet_sta(sind:eind);
    
    % CDEC FNF for FOL
    cdec_fnf_sta = load(['G:\CA-LSTM\CDEC\FNF_mm_QC\FNF_',cdec_id{i},'_mm.txt']); % [yr mon day flow(mm)]
    sind = find(cdec_fnf_sta(:,1)==val_startdate(1)&cdec_fnf_sta(:,2)==val_startdate(2)&cdec_fnf_sta(:,3)==val_startdate(3));
    eind = find(cdec_fnf_sta(:,1)==val_enddate(1)&cdec_fnf_sta(:,2)==val_enddate(2)&cdec_fnf_sta(:,3)==val_enddate(3));
    cdec_fnf_sta_val = cdec_fnf_sta(sind:eind,4);
    
    % NSE: any negative FLOW_cal_obs values (either missing or FNF calculation) are excluded
    nse = 1 - (mean((simflow_sta_val(cdec_fnf_sta_val>=0)-cdec_fnf_sta_val(cdec_fnf_sta_val>=0)).^2)/var(cdec_fnf_sta_val(cdec_fnf_sta_val>=0),1));
    
    % KGE
    mean_ratio = mean(simflow_sta_val(cdec_fnf_sta_val>=0))/mean(cdec_fnf_sta_val(cdec_fnf_sta_val>=0));
    std_ratio = std(simflow_sta_val(cdec_fnf_sta_val>=0),1)/std(cdec_fnf_sta_val(cdec_fnf_sta_val>=0),1);
    lin_corr = corr(simflow_sta_val(cdec_fnf_sta_val>=0),cdec_fnf_sta_val(cdec_fnf_sta_val>=0));
    kge = 1-((1-mean_ratio)^2+(1-std_ratio)^2+(1-lin_corr)^2)^0.5;
    
    
    kge_nse_val(i,1) = kge;
    kge_nse_val(i,2) = nse;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savedata_sta_cal = [cal_datemat_cal(:,1:3),simflow_sta_cal,cdec_fnf_sta_cal];
    savedata_sta_val = [cal_datemat_val(:,1:3),simflow_sta_val,cdec_fnf_sta_val];
    savedata_sta = [savedata_sta_cal;savedata_sta_val];
%     fid = fopen(['G:\CA-LSTM\HYMOD\calibration\',trial_sel,'\simflow_watershed\simflow_hymod_',cdec_id{i},'.txt'],'w');
%     fprintf(fid,'%4d %2d %2d %14.8f %14.8f\n',savedata_sta');
%     fclose(fid);
    
    %aet_sta = [savedata_sta(:,1:3), [simaet_sta_cal;simaet_sta_val]];
    aet_sta = [sim_datemat(:,1:3), simaet_sta];
%     fid = fopen(['G:\CA-LSTM\HYMOD\calibration\',trial_sel,'\simaet_watershed\simaet_hymod_',cdec_id{i},'.txt'],'w');
%     fprintf(fid,'%4d %2d %2d %14.8f\n',aet_sta');
%     fclose(fid);
	
	swe_sta = [sim_datemat(:,1:3), simswe_sta];
    fid = fopen(['G:\CA-LSTM\HYMOD\calibration\',trial_sel,'\simswe_watershed\simswe_hymod_',cdec_id{i},'.txt'],'w');
    fprintf(fid,'%4d %2d %2d %14.8f\n',swe_sta');
    fclose(fid);
    
end

% save(['G:\CA-LSTM\HYMOD\calibration\',trial_sel,'\kge_nse_cal.txt'],'kge_nse_cal','-ascii');
% save(['G:\CA-LSTM\HYMOD\calibration\',trial_sel,'\kge_nse_val.txt'],'kge_nse_val','-ascii');

end


