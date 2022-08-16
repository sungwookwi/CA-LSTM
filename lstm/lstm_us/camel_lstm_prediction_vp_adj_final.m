
clear
clc

%tic;
dT_all = [0,1,2,3,4];

basin_list = importdata('./basin_list.txt','\n',531);
gauge_info = xlsread('./gauge_information_531.xlsx');

% SPLIT DATA INTO TRAIN (1999-2008) & TEST (1989-1999)
eyr_train = 2008;
syr_test = 1988;

% MAURER CLIMATE for each basin
datemat_maurer = datevec(datenum([1980 1 1]):datenum([2008 12 31]));
camels_tmin = nan(datenum([eyr_train,12,31])-datenum([syr_test,1,1])+1,length(basin_list));
camels_vp = nan(datenum([eyr_train,12,31])-datenum([syr_test,1,1])+1,length(basin_list));

f_power_all = nan(length(basin_list),3);
gof_power_all = nan(length(basin_list),1);
for i = 1:length(basin_list)
    
    basin_id = basin_list{i};
    if gauge_info(gauge_info(:,2)==str2double(basin_id),1)<10
        basin_huc = ['0',num2str(gauge_info(gauge_info(:,2)==str2double(basin_id),1))];
    else
        basin_huc = num2str(gauge_info(gauge_info(:,2)==str2double(basin_id),1));
    end
    
    basin_camels = importdata(['./maurer_extended/',basin_huc,'/',basin_id,'_lump_maurer_forcing_leap.txt'],'\t',4);

    % basin_camels.colheaders : Year Mnth Day Hr Dayl(s) PRCP(mm/day)
    % SRAD(W/m2) SWE(mm) Tmax(C) Tmin(C) Vp(Pa)
    basin_tmin = basin_camels.data(:,10)+273.15;
    basin_vp = basin_camels.data(:,11);   
    
    sind_all = find(datemat_maurer(:,1)==syr_test & datemat_maurer(:,2)==1 & datemat_maurer(:,3)==1);
    eind_all = find(datemat_maurer(:,1)==eyr_train & datemat_maurer(:,2)==12 & datemat_maurer(:,3)==31);
    basin_tmin_all = basin_tmin(sind_all:eind_all);
    basin_vp_all = basin_vp(sind_all:eind_all);
    
	[f_power,gof_power]=fit(basin_tmin_all,basin_vp_all,'power2');
    
    if gof_power.rsquare >= .8    
        f_power_all(i,:) = [f_power.a,f_power.b,f_power.c];
    else
        camels_tmin(:,i) = basin_tmin_all;
        camels_vp(:,i) = basin_vp_all;
    end
    gof_power_all(i) = gof_power.rsquare;   
 
end

camels_tmin_remain = camels_tmin(:);
camels_vp_remain = camels_vp(:);
camels_tmin_remain(isnan(camels_tmin_remain))=[];
camels_vp_remain(isnan(camels_vp_remain))=[];

[f_power_remain,~]=fit(camels_tmin_remain,camels_vp_remain,'power2');

nan_ind = find(isnan(f_power_all(:,1)));
for i = 1:length(nan_ind)
    f_power_all(nan_ind(i),:) = [f_power_remain.a,f_power_remain.b,f_power_remain.c];
end


for ttttt = 7%8:10

trial = ttttt;


for ttt = 4:5


dT = dT_all(ttt);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       LOAD CAMELS BASIN DATA (531 basins)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;

%----------------------------- MAURER CLIMATE ------------------------------------
basin_list = importdata('./basin_list.txt','\n',531);
gauge_info = xlsread('./gauge_information_531.xlsx');
basin_list_num = str2double(basin_list);

% SPLIT DATA INTO TRAIN (1999-2008) & TEST (1989-1999)
syr_train = 1998;
eyr_train = 2008;
syr_test = 1988;
eyr_test = 1999;

% MAURER CLIMATE for each basin
datemat_maurer = datevec(datenum([1980 1 1]):datenum([2008 12 31]));
datemat_maurer_train = datevec(datenum([syr_train 1 1]):datenum([eyr_train 12 31]));
datemat_maurer_test = datevec(datenum([syr_test 1 1]):datenum([eyr_test 12 31]));
basin_area = nan(length(basin_list),1);
for i = 1:length(basin_list)
    
    basin_id = basin_list{i};
    if gauge_info(gauge_info(:,2)==str2double(basin_id),1)<10
        basin_huc = ['0',num2str(gauge_info(gauge_info(:,2)==str2double(basin_id),1))];
    else
        basin_huc = num2str(gauge_info(gauge_info(:,2)==str2double(basin_id),1));
    end
    
    basin_camels = importdata(['./maurer_extended/',basin_huc,'/',basin_id,'_lump_maurer_forcing_leap.txt'],'\t',4);
    basin_area(i) = str2double(basin_camels.textdata{3});
    % basin_camels.colheaders : Year Mnth Day Hr Dayl(s) PRCP(mm/day)
    % SRAD(W/m2) SWE(mm) Tmax(C) Tmin(C) Vp(Pa)
    basin_pr = basin_camels.data(:,6);
    basin_tmax = basin_camels.data(:,9)+dT;
    basin_tmin = basin_camels.data(:,10)+dT;
    basin_srad = basin_camels.data(:,7);
    %basin_vp = basin_camels.data(:,11);

basin_tmin_K = basin_tmin +273.15;
    basin_tmin_K_org = basin_camels.data(:,10)+273.15;
    basin_vp_org = basin_camels.data(:,11);
    
    vp_delta = (f_power_all(i,1)*basin_tmin_K.^f_power_all(i,2)+f_power_all(i,3))-(f_power_all(i,1)*basin_tmin_K_org.^f_power_all(i,2)+f_power_all(i,3));
    basin_vp = basin_vp_org+vp_delta;
	
    
    sind_train = find(datemat_maurer(:,1)==syr_train & datemat_maurer(:,2)==1 & datemat_maurer(:,3)==1);
    eind_train = find(datemat_maurer(:,1)==eyr_train & datemat_maurer(:,2)==12 & datemat_maurer(:,3)==31);
    sind_test = find(datemat_maurer(:,1)==syr_test & datemat_maurer(:,2)==1 & datemat_maurer(:,3)==1);
    eind_test = find(datemat_maurer(:,1)==eyr_test & datemat_maurer(:,2)==12 & datemat_maurer(:,3)==31);         
    
    basin_pr_train = basin_pr(sind_train:eind_train);
    basin_tmax_train = basin_tmax(sind_train:eind_train);
    basin_tmin_train = basin_tmin(sind_train:eind_train);
    basin_srad_train = basin_srad(sind_train:eind_train);
    basin_vp_train = basin_vp(sind_train:eind_train);
    
    basin_pr_test = basin_pr(sind_test:eind_test);
    basin_tmax_test = basin_tmax(sind_test:eind_test);
    basin_tmin_test = basin_tmin(sind_test:eind_test);
    basin_srad_test = basin_srad(sind_test:eind_test);
    basin_vp_test = basin_vp(sind_test:eind_test);
    
    
    eval(['pr_train_',basin_id,' = basin_pr_train;'])
    eval(['tmax_train_',basin_id,' = basin_tmax_train;'])
    eval(['tmin_train_',basin_id,' = basin_tmin_train;'])
    eval(['srad_train_',basin_id,' = basin_srad_train;'])
    eval(['vp_train_',basin_id,' = basin_vp_train;'])
    
    eval(['pr_test_',basin_id,' = basin_pr_test;'])
    eval(['tmax_test_',basin_id,' = basin_tmax_test;'])
    eval(['tmin_test_',basin_id,' = basin_tmin_test;'])
    eval(['srad_test_',basin_id,' = basin_srad_test;'])
    eval(['vp_test_',basin_id,' = basin_vp_test;'])
    
end


[dyn_mean_std,dyn_name] = xlsread(['./train_data_scaler/train_data_scaler_t',num2str(trial),'.xls'],'dyn');
name_ind = find(ismember(dyn_name,'PRCP'));
basin_pr_train_all_mean = dyn_mean_std(name_ind,1);
basin_pr_train_all_std = dyn_mean_std(name_ind,2);
name_ind = find(ismember(dyn_name,'Tmax'));
basin_tmax_train_all_mean = dyn_mean_std(name_ind,1);
basin_tmax_train_all_std = dyn_mean_std(name_ind,2);
name_ind = find(ismember(dyn_name,'Tmin'));
basin_tmin_train_all_mean = dyn_mean_std(name_ind,1);
basin_tmin_train_all_std = dyn_mean_std(name_ind,2);
name_ind = find(ismember(dyn_name,'SRAD'));
basin_srad_train_all_mean = dyn_mean_std(name_ind,1);
basin_srad_train_all_std = dyn_mean_std(name_ind,2);
name_ind = find(ismember(dyn_name,'Vp'));
basin_vp_train_all_mean = dyn_mean_std(name_ind,1);
basin_vp_train_all_std = dyn_mean_std(name_ind,2);


% Standardizing MAURER CLIMATE for each basin
% e.g., pr_standard_train_01022500; pr_standard_test_01022500; tmax_standard_train_01022500; tmax_standard_test_01022500
for i = 1:length(basin_list)
    eval(['basin_pr_train = pr_train_',basin_list{i},';'])
    eval(['pr_standard_train_',basin_list{i},' = (basin_pr_train-basin_pr_train_all_mean)/basin_pr_train_all_std;'])
    eval(['clear(''pr_train_',basin_list{i},''')'])
    eval(['basin_pr_test = pr_test_',basin_list{i},';'])
    eval(['pr_standard_test_',basin_list{i},' = (basin_pr_test-basin_pr_train_all_mean)/basin_pr_train_all_std;'])
    eval(['clear(''pr_test_',basin_list{i},''')'])
    
    eval(['basin_tmax_train = tmax_train_',basin_list{i},';'])
    eval(['tmax_standard_train_',basin_list{i},' = (basin_tmax_train-basin_tmax_train_all_mean)/basin_tmax_train_all_std;'])
    eval(['clear(''tmax_train_',basin_list{i},''')'])
    eval(['basin_tmax_test = tmax_test_',basin_list{i},';'])
    eval(['tmax_standard_test_',basin_list{i},' = (basin_tmax_test-basin_tmax_train_all_mean)/basin_tmax_train_all_std;'])
    eval(['clear(''tmax_test_',basin_list{i},''')'])
    
    eval(['basin_tmin_train = tmin_train_',basin_list{i},';'])
    eval(['tmin_standard_train_',basin_list{i},' = (basin_tmin_train-basin_tmin_train_all_mean)/basin_tmin_train_all_std;'])
    eval(['clear(''tmin_train_',basin_list{i},''')'])
    eval(['basin_tmin_test = tmin_test_',basin_list{i},';'])
    eval(['tmin_standard_test_',basin_list{i},' = (basin_tmin_test-basin_tmin_train_all_mean)/basin_tmin_train_all_std;'])
    eval(['clear(''tmin_test_',basin_list{i},''')'])
    
    eval(['basin_srad_train = srad_train_',basin_list{i},';'])
    eval(['srad_standard_train_',basin_list{i},' = (basin_srad_train-basin_srad_train_all_mean)/basin_srad_train_all_std;'])
    eval(['clear(''srad_train_',basin_list{i},''')'])
    eval(['basin_srad_test = srad_test_',basin_list{i},';'])
    eval(['srad_standard_test_',basin_list{i},' = (basin_srad_test-basin_srad_train_all_mean)/basin_srad_train_all_std;'])
    eval(['clear(''srad_test_',basin_list{i},''')'])
    
    eval(['basin_vp_train = vp_train_',basin_list{i},';'])
    eval(['vp_standard_train_',basin_list{i},' = (basin_vp_train-basin_vp_train_all_mean)/basin_vp_train_all_std;'])
    eval(['clear(''vp_train_',basin_list{i},''')'])
    eval(['basin_vp_test = vp_test_',basin_list{i},';'])
    eval(['vp_standard_test_',basin_list{i},' = (basin_vp_test-basin_vp_train_all_mean)/basin_vp_train_all_std;'])
    eval(['clear(''vp_test_',basin_list{i},''')'])
    
end


%------------------------ STATIC ATTRIBUTE --------------------------------
[static_mean_std,static_name] = xlsread(['./train_data_scaler/train_data_scaler_t',num2str(trial),'.xls'],'camel_attr');

% camels_clim
[camels_clim,var_clim] = xlsread('./camels_attributes_v2.0/camels_attr_all.xlsx','camels_clim');
basin_attr_p_mean = nan(length(basin_list),1);
basin_attr_pet_mean = nan(length(basin_list),1);
basin_attr_aridity = nan(length(basin_list),1);
basin_attr_frac_snow = nan(length(basin_list),1);
basin_attr_high_prec_freq = nan(length(basin_list),1);
basin_attr_high_prec_dur = nan(length(basin_list),1);
basin_attr_low_prec_freq = nan(length(basin_list),1);
basin_attr_low_prec_dur = nan(length(basin_list),1);
for i = 1:length(basin_list)
    % Precipitation mean
    basin_attr_p_mean(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'p_mean'));
    % PET mean
    basin_attr_pet_mean(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'pet_mean'));
    % Aridity index
    basin_attr_aridity(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'aridity'));
    % Snow fraction
    basin_attr_frac_snow(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'frac_snow'));
    % High precip frequency
    basin_attr_high_prec_freq(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'high_prec_freq'));
    % High precip duration
    basin_attr_high_prec_dur(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'high_prec_dur'));
    % Low precip frequency
    basin_attr_low_prec_freq(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'low_prec_freq'));
    % Low precip duration
    basin_attr_low_prec_dur(i) = camels_clim(camels_clim(:,1)==str2double(basin_list{i}),ismember(var_clim(1,:),'low_prec_dur'));
end

name_ind = find(ismember(static_name,'p_mean'));
basin_attr_p_mean_standard = (basin_attr_p_mean-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'pet_mean'));
basin_attr_pet_mean_standard = (basin_attr_pet_mean-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'aridity'));
basin_attr_aridity_standard = (basin_attr_aridity-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'frac_snow'));
basin_attr_frac_snow_standard = (basin_attr_frac_snow-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'high_prec_freq'));
basin_attr_high_prec_freq_standard = (basin_attr_high_prec_freq-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'high_prec_dur'));
basin_attr_high_prec_dur_standard = (basin_attr_high_prec_dur-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'low_prec_freq'));
basin_attr_low_prec_freq_standard = (basin_attr_low_prec_freq-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'low_prec_dur'));
basin_attr_low_prec_dur_standard = (basin_attr_low_prec_dur-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);

% camels_topo
[camels_topo,var_topo] = xlsread('./camels_attributes_v2.0/camels_attr_all.xlsx','camels_topo');
basin_attr_elev_mean = nan(length(basin_list),1);
basin_attr_slope_mean = nan(length(basin_list),1);
basin_attr_area_gages2 = nan(length(basin_list),1);
for i = 1:length(basin_list)
    % Elevation
    basin_attr_elev_mean(i) = camels_topo(camels_topo(:,1)==str2double(basin_list{i}),ismember(var_topo(1,:),'elev_mean'));
    % Slope
    basin_attr_slope_mean(i) = camels_topo(camels_topo(:,1)==str2double(basin_list{i}),ismember(var_topo(1,:),'slope_mean'));
    % Area
    basin_attr_area_gages2(i) = camels_topo(camels_topo(:,1)==str2double(basin_list{i}),ismember(var_topo(1,:),'area_gages2'));
end
name_ind = find(ismember(static_name,'elev_mean'));
basin_attr_elev_mean_standard = (basin_attr_elev_mean-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'slope_mean'));
basin_attr_slope_mean_standard = (basin_attr_slope_mean-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'area_gages2'));
basin_attr_area_gages2_standard = (basin_attr_area_gages2-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);

% camels_vege
[camels_vege,var_vege] = xlsread('./camels_attributes_v2.0/camels_attr_all.xlsx','camels_vege');
basin_attr_frac_forest = nan(length(basin_list),1);
basin_attr_lai_max = nan(length(basin_list),1);
basin_attr_lai_diff = nan(length(basin_list),1);
basin_attr_gvf_max = nan(length(basin_list),1);
basin_attr_gvf_diff = nan(length(basin_list),1);
for i = 1:length(basin_list)
    % Forest fraction
    basin_attr_frac_forest(i) = camels_vege(camels_vege(:,1)==str2double(basin_list{i}),ismember(var_vege(1,:),'frac_forest'));
    % LAI max
    basin_attr_lai_max(i) = camels_vege(camels_vege(:,1)==str2double(basin_list{i}),ismember(var_vege(1,:),'lai_max'));
    % LAI difference
    basin_attr_lai_diff(i) = camels_vege(camels_vege(:,1)==str2double(basin_list{i}),ismember(var_vege(1,:),'lai_diff'));
    % GVF max
    basin_attr_gvf_max(i) = camels_vege(camels_vege(:,1)==str2double(basin_list{i}),ismember(var_vege(1,:),'gvf_max'));
    % GVF difference
    basin_attr_gvf_diff(i) = camels_vege(camels_vege(:,1)==str2double(basin_list{i}),ismember(var_vege(1,:),'gvf_diff'));
end
name_ind = find(ismember(static_name,'frac_forest'));
basin_attr_frac_forest_standard = (basin_attr_frac_forest-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'lai_max'));
basin_attr_lai_max_standard = (basin_attr_lai_max-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'lai_diff'));
basin_attr_lai_diff_standard = (basin_attr_lai_diff-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'gvf_max'));
basin_attr_gvf_max_standard = (basin_attr_gvf_max-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'gvf_diff'));
basin_attr_gvf_diff_standard = (basin_attr_gvf_diff-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);

% camels_soil
[camels_soil,var_soil] = xlsread('./camels_attributes_v2.0/camels_attr_all.xlsx','camels_soil');
basin_attr_soil_depth_pelletier = nan(length(basin_list),1);
basin_attr_soil_depth_statsgo = nan(length(basin_list),1);
basin_attr_soil_porosity = nan(length(basin_list),1);
basin_attr_soil_conductivity = nan(length(basin_list),1);
basin_attr_max_water_content = nan(length(basin_list),1);
basin_attr_sand_frac = nan(length(basin_list),1);
basin_attr_silt_frac = nan(length(basin_list),1);
basin_attr_clay_frac = nan(length(basin_list),1);
for i = 1:length(basin_list)
    % Soil depth (Pelletier)
    basin_attr_soil_depth_pelletier(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'soil_depth_pelletier'));
    % Soil depth (STATSGO)
    basin_attr_soil_depth_statsgo(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'soil_depth_statsgo'));
    % Soil porosity
    basin_attr_soil_porosity(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'soil_porosity'));
    % Soil conductivity
    basin_attr_soil_conductivity(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'soil_conductivity'));
    % Max water content
    basin_attr_max_water_content(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'max_water_content'));
    % Sand fraction
    basin_attr_sand_frac(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'sand_frac'));
    % Silt fraction
    basin_attr_silt_frac(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'silt_frac'));
    % Clay fraction
    basin_attr_clay_frac(i) = camels_soil(camels_soil(:,1)==str2double(basin_list{i}),ismember(var_soil(1,:),'clay_frac'));
end
name_ind = find(ismember(static_name,'soil_depth_pelletier'));
basin_attr_soil_depth_pelletier_standard = (basin_attr_soil_depth_pelletier-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'soil_depth_statsgo'));
basin_attr_soil_depth_statsgo_standard = (basin_attr_soil_depth_statsgo-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'soil_porosity'));
basin_attr_soil_porosity_standard = (basin_attr_soil_porosity-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'soil_conductivity'));
basin_attr_soil_conductivity_standard = (basin_attr_soil_conductivity-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'max_water_content'));
basin_attr_max_water_content_standard = (basin_attr_max_water_content-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'sand_frac'));
basin_attr_sand_frac_standard = (basin_attr_sand_frac-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'silt_frac'));
basin_attr_silt_frac_standard = (basin_attr_silt_frac-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'clay_frac'));
basin_attr_clay_frac_standard = (basin_attr_clay_frac-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);

% camels_geol
[camels_geol,var_geol] = xlsread('./camels_attributes_v2.0/camels_attr_all.xlsx','camels_geol');
basin_attr_carbonate_rocks_frac = nan(length(basin_list),1);
basin_attr_geol_permeability = nan(length(basin_list),1);
for i = 1:length(basin_list)
    % Carbonate rocks fraction
    basin_attr_carbonate_rocks_frac(i) = camels_geol(camels_geol(:,1)==str2double(basin_list{i}),ismember(var_geol(1,:),'carbonate_rocks_frac'));
    % Geological permeability
    basin_attr_geol_permeability(i) = camels_geol(camels_geol(:,1)==str2double(basin_list{i}),ismember(var_geol(1,:),'geol_permeability'));
end
name_ind = find(ismember(static_name,'carbonate_rocks_frac'));
basin_attr_carbonate_rocks_frac_standard = (basin_attr_carbonate_rocks_frac-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);
name_ind = find(ismember(static_name,'geol_permeability'));
basin_attr_geol_permeability_standard = (basin_attr_geol_permeability-static_mean_std(name_ind,1))/static_mean_std(name_ind,2);


%------------------------ STREAMFLOW --------------------------------

% CAMELS STREAMFLOW for each basin
datemat_train = datevec(datenum([1999 10 1]):datenum([2008 9 30]));
datemat_test = datevec(datenum([1989 10 1]):datenum([1999 9 30]));
basin_all_id_std_size = nan(length(basin_list),3);
basin_all_id_std_size_test = nan(length(basin_list),3);
for i = 1:length(basin_list)   
    
    basin_id = basin_list{i};
    basin_data = importdata(['./usgs_streamflow_reformat/',basin_id,'_streamflow_qc.txt'],' ');    
        
    sind_train = find(basin_data(:,1)==datemat_train(1,1) & basin_data(:,2)==datemat_train(1,2) & basin_data(:,3)==datemat_train(1,3));
    eind_train = find(basin_data(:,1)==datemat_train(end,1) & basin_data(:,2)==datemat_train(end,2) & basin_data(:,3)==datemat_train(end,3));
    sind_test = find(basin_data(:,1)==datemat_test(1,1) & basin_data(:,2)==datemat_test(1,2) & basin_data(:,3)==datemat_test(1,3));
    eind_test = find(basin_data(:,1)==datemat_test(end,1) & basin_data(:,2)==datemat_test(end,2) & basin_data(:,3)==datemat_test(end,3));  
    if isempty(sind_train)
        datenum_first = datenum(basin_data(1,1:3));
        datemat_miss = datevec(datenum(datemat_train(1,:)):datenum_first-1);
        basin_q_train = [[datemat_miss(:,1:3),repmat(-999,size(datemat_miss,1),1)];basin_data(1:eind_train,:)];
    else
        basin_q_train = basin_data(sind_train:eind_train,:);
    end
    if isempty(eind_test)
        datenum_last = datenum(basin_data(end,1:3));
        datemat_miss = datevec(datenum_last+1:datenum(datemat_test(end,:)));
        basin_q_test = [basin_data(sind_test:end,:);[datemat_miss(:,1:3),repmat(-999,size(datemat_miss,1),1)]];
    else
        basin_q_test = basin_data(sind_test:eind_test,:);
    end
    

    basin_q_train(basin_q_train(:,4)<0,4) = nan;
    basin_q_train(:,4) = basin_q_train(:,4)*0.3048^3*3600*24/basin_area(i)*1000;
    basin_q_test(basin_q_test(:,4)<0,4) = nan;
    basin_q_test(:,4) = basin_q_test(:,4)*0.3048^3*3600*24/basin_area(i)*1000;

    
    basin_all_id_std_size(i,:) = [str2double(basin_id),nanstd(basin_q_train(:,4),1),sum(~isnan(basin_q_train(:,4)))];
    basin_all_id_std_size_test(i,:) = [str2double(basin_id),nanstd(basin_q_test(:,4),1),sum(~isnan(basin_q_test(:,4)))];
    
    eval(['q_train_',basin_id,' = basin_q_train;'])
    eval(['q_test_',basin_id,' = basin_q_test;'])
    
end

q_mean_std = xlsread(['./train_data_scaler/train_data_scaler_t',num2str(trial),'.xls'],'target');


tic;
for ii = 1:size(basin_all_id_std_size,1)
    
    basin_id = basin_list{ii};
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                    LSTM INPUT FOR TRAINING & TEST
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    num_feature_dyn = 5;
    seq_length = 365;
    
%     % TRAIN
%     tot_sample_num = sum(basin_all_id_std_size(ii,end));
%     XTrain_basin_id = nan(tot_sample_num,1);
%     XTrain = cell(tot_sample_num,1);
%     YTrain = cell(tot_sample_num,1);
%     
%     
%     eval(['q_basin = q_train_',basin_id,';'])
%     
%     eval(['pr_basin = pr_standard_train_',basin_id,';'])
%     eval(['tmax_basin = tmax_standard_train_',basin_id,';'])
%     eval(['tmin_basin = tmin_standard_train_',basin_id,';'])
%     eval(['srad_basin = srad_standard_train_',basin_id,';'])
%     eval(['vp_basin = vp_standard_train_',basin_id,';'])
%     
%     n=1;
%     for j = 1:size(q_basin,1)
%         
%         if ~isnan(q_basin(j,end))
%             
%             meteo_ind = find(datemat_maurer_train(:,1)==q_basin(j,1)&datemat_maurer_train(:,2)==q_basin(j,2)&datemat_maurer_train(:,3)==q_basin(j,3));
%             x_seq = nan(num_feature_dyn,seq_length);
%             nn = 1;
%             for tt = meteo_ind-seq_length+1:meteo_ind
%                 x_seq(:,nn) = [pr_basin(tt);srad_basin(tt);tmax_basin(tt)+dT;tmin_basin(tt)+dT;vp_basin(tt);...
%                     ];
%                 nn=nn+1;
%             end
%             XTrain{n} = x_seq;
%             YTrain{n} = q_basin(j,end);
%             XTrain_basin_id(n) = str2double(basin_id);
%             n=n+1;
%         end
%     end
    
    
    % TEST
    tot_sample_num_test = sum(basin_all_id_std_size_test(ii,end));
    XTest_basin_id = nan(tot_sample_num_test,1);
    XTest = cell(tot_sample_num_test,1);
    YTest = cell(tot_sample_num_test,1);
    n=1;
    
    eval(['q_basin = q_test_',basin_id,';'])
    
    eval(['pr_basin = pr_standard_test_',basin_id,';'])
    eval(['tmax_basin = tmax_standard_test_',basin_id,';'])
    eval(['tmin_basin = tmin_standard_test_',basin_id,';'])
    eval(['srad_basin = srad_standard_test_',basin_id,';'])
    eval(['vp_basin = vp_standard_test_',basin_id,';'])
    
    for j = 1:size(q_basin,1)
        
        if ~isnan(q_basin(j,end))
            
            meteo_ind = find(datemat_maurer_test(:,1)==q_basin(j,1)&datemat_maurer_test(:,2)==q_basin(j,2)&datemat_maurer_test(:,3)==q_basin(j,3));
            x_seq = nan(num_feature_dyn,seq_length);
            nn = 1;
            for tt = meteo_ind-seq_length+1:meteo_ind
                x_seq(:,nn) = [pr_basin(tt);srad_basin(tt);tmax_basin(tt);tmin_basin(tt);vp_basin(tt);...
                    ];
                nn=nn+1;
            end
            XTest{n} = x_seq;
            YTest{n} = q_basin(j,end);
            XTest_basin_id(n) = str2double(basin_id);
            n=n+1;
        end
    end
        
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                             LSTM PREDICTION
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %------------------- Hyperparameter not being tuned -----------------------
    % Number of input features
    num_feature_static = 26;
    num_input_lyr_unit = num_feature_dyn + num_feature_static;
    
    % Number of hidden layer units
    num_hiddn_lyr_unit_1 = 256;
    
    % Row indies for each gate weights
    iInd = (1:num_hiddn_lyr_unit_1*1);
    fInd = (num_hiddn_lyr_unit_1*1+1:num_hiddn_lyr_unit_1*2);
    gInd = (num_hiddn_lyr_unit_1*2+1:num_hiddn_lyr_unit_1*3);
    oInd = (num_hiddn_lyr_unit_1*3+1:num_hiddn_lyr_unit_1*4);
    ifoInd = [iInd fInd oInd];
    
    
    %------------------- LSTM Weight Loading ---------------------------
    % Input Weight
    U = load(['./weights/lstm_weight_ih_t',num2str(trial),'.txt']);
    % Hidden Weight
    W = load(['./weights/lstm_weight_hh_t',num2str(trial),'.txt']);
    % Output Weight
    V = load(['./weights/lstm_weight_ho_t',num2str(trial),'.txt']);
    % Bias
    b_ih = load(['./weights/lstm_bias_ih_t',num2str(trial),'.txt']);
    b_hh = load(['./weights/lstm_bias_hh_t',num2str(trial),'.txt']);
    bv = load(['./weights/lstm_bias_ho_t',num2str(trial),'.txt']);
    
    U = gpuArray(U);
    W = gpuArray(W);
    V = gpuArray(V);
    b_ih = gpuArray(b_ih);
    b_hh = gpuArray(b_hh);
    bv = gpuArray(bv);
    
    
    %--------------------------- LSTM Prediction ------------------------------

    attr_static_basin = [...
        basin_attr_area_gages2_standard(ii);...
        basin_attr_aridity_standard(ii);...
        basin_attr_carbonate_rocks_frac_standard(ii);...
        basin_attr_clay_frac_standard(ii);...
        basin_attr_elev_mean_standard(ii);...
        basin_attr_frac_forest_standard(ii);...
        basin_attr_frac_snow_standard(ii);...
        basin_attr_geol_permeability_standard(ii);... 
        basin_attr_gvf_diff_standard(ii);...
        basin_attr_gvf_max_standard(ii);...
        basin_attr_high_prec_dur_standard(ii);...
        basin_attr_high_prec_freq_standard(ii);...
        basin_attr_lai_diff_standard(ii);...
        basin_attr_lai_max_standard(ii);...
        basin_attr_low_prec_dur_standard(ii);...
        basin_attr_low_prec_freq_standard(ii);...
        basin_attr_max_water_content_standard(ii);...
        basin_attr_p_mean_standard(ii);...
        basin_attr_pet_mean_standard(ii);...       
        basin_attr_sand_frac_standard(ii);...
        basin_attr_silt_frac_standard(ii);...
        basin_attr_slope_mean_standard(ii);... 
        basin_attr_soil_conductivity_standard(ii);...
        basin_attr_soil_depth_pelletier_standard(ii);...
        basin_attr_soil_depth_statsgo_standard(ii);...
        basin_attr_soil_porosity_standard(ii);...   
        ];
    %%%%%%%%%%%%%%%%%%%%%%%%% TRAIN PERIOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     minibatch_size = size(XTrain,1);
%     
%     X_static = reshape(repmat(attr_static_basin,minibatch_size,seq_length),num_feature_static,minibatch_size,seq_length);
%     X_dynamic = permute(reshape([XTrain{:}],num_feature_dyn,seq_length,minibatch_size),[1,3,2]);
%     
%     X_all = zeros(num_input_lyr_unit,minibatch_size,seq_length,'single');
%     X_all(1:num_feature_dyn,:,:) = X_dynamic;
%     X_all(num_feature_dyn+1:num_input_lyr_unit,:,:) = X_static;
%    
%     
%     tic;
%     % FORWARD
%     X =X_all;
%     H_all = zeros(num_hiddn_lyr_unit_1,minibatch_size,seq_length,'single');
%     C_all = zeros(num_hiddn_lyr_unit_1,minibatch_size,seq_length,'single');
%     H = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single');
%     C = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single');
%     
%     for tt = 1:seq_length
%         % Gate operation
%         G = U*X(:,:,tt)+b_ih + W*H + b_hh;
%         G(ifoInd,:) = 1./(1+exp(-G(ifoInd,:)));
%         G(gInd,:)=tanh(G(gInd,:));
%         
%         %Cell state update
%         C = G(gInd,:).*G(iInd,:)+G(fInd,:).*C;
%         C_all(:,:,tt) = C;
%         
%         % LSTM layer output
%         H = tanh(C).*G(oInd, :);
%         H_all(:,:,tt) = H;
%         
%     end
%     toc;
%     
%     % Output layer output
%     Yhat_train = (V*H + bv)*q_mean_std(2)+q_mean_std(1);
%     Yhat_train(1,Yhat_train(1,:)<0) = 0;
%     Yobs_train = single(reshape([YTrain{:}],1,minibatch_size));
%     
%     Yobs_valid = Yobs_train(Yobs_train>=0);
%     Yhat_valid = Yhat_train(Yobs_train>=0);
%     
%     nse_train = 1-mean((Yhat_valid-Yobs_valid).^2)/var(Yobs_valid,1);
    %nse_train_all(j) = double(nse_train);
    
    
%     figure
%     plot(Yhat_train,'k'); hold on;
%     plot(Yobs_train,'b');
    
    %     figure
    %     plot(Yobs,Yhat,'m*')
    
    
    %tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%% TEST PERIOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minibatch_size = size(XTest,1);
    
    X_static = reshape(repmat(attr_static_basin,minibatch_size,seq_length),num_feature_static,minibatch_size,seq_length);
    X_dynamic = permute(reshape([XTest{:}],num_feature_dyn,seq_length,minibatch_size),[1,3,2]);
    
    X_all = zeros(num_input_lyr_unit,minibatch_size,seq_length,'single','gpuArray');
    X_all(1:num_feature_dyn,:,:) = X_dynamic;
    X_all(num_feature_dyn+1:num_input_lyr_unit,:,:) = X_static;

    
    
    % FORWARD
    X =X_all;
    %H_all = zeros(num_hiddn_lyr_unit_1,minibatch_size,seq_length,'single','gpuArray');
    %C_all = zeros(num_hiddn_lyr_unit_1,minibatch_size,seq_length,'single','gpuArray');
    H = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');
    C = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');
    
    for tt = 1:seq_length
        % Gate operation
        G = U*X(:,:,tt)+b_ih + W*H + b_hh;
        G(ifoInd,:) = 1./(1+exp(-G(ifoInd,:)));
        G(gInd,:)=tanh(G(gInd,:));
        
        %Cell state update
        C = G(gInd,:).*G(iInd,:)+G(fInd,:).*C;
        %C_all(:,:,tt) = C;
        
        % LSTM layer output
        H = tanh(C).*G(oInd, :);
        %H_all(:,:,tt) = H;
        
    end
    %toc;
    
    % Output layer output
    Yhat_test = (V*H + bv)*q_mean_std(2)+q_mean_std(1);
    Yhat_test(1,Yhat_test(1,:)<0) = 0;
    Yobs_test = single(reshape([YTest{:}],1,minibatch_size));
    
    Yobs_valid = Yobs_test(Yobs_test>=0);
    Yhat_valid = Yhat_test(Yobs_test>=0);
    
    nse_test = 1-mean((Yhat_valid-Yobs_valid).^2)/var(Yobs_valid,1);
    %nse_test_all(j) = double(nse_test);
    
%     figure
%     plot(Yhat_test,'k'); hold on;
%     plot(Yobs_test,'b');
    
    %     figure
    %     plot(Yobs,Yhat,'m*')
    
    Yhat_test_save = gather(Yhat_test);

    savedata = double([Yhat_test_save',Yobs_test']);
    save(['./Qsim/trial_',num2str(trial),'/dT',num2str(dT),'/SimQ_ObsQ_',basin_list{ii},'.txt'],'savedata','-ascii')
    
e_t = toc;
    disp(['trial ',num2str(trial),' dT ',num2str(dT),' basin# ',num2str(ii),'  Elapsed Time:  ',num2str(e_t),' sec'])
end
%[nse_calib_train_all,nse_calib_test_all]
toc;


end
end