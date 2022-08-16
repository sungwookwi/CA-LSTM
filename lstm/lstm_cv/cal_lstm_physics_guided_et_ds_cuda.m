% Long Short Term Memory Network for Predictions in California
% Written by Sungwook Wi

clear
clc

% w_q = [0.1, 0.3, 0.5, 0.7, 0.9];
% w_et  = [0.9, 0.7, 0.5, 0.3, 0.1];
% alpha_str = {'01','03','05','07','09'};

w_q = [0.1, 0.3, 0.5, 0.7, 0.9];
w_et  = [0.9, 0.7, 0.5, 0.3, 0.1];
alpha_str = {'01','03','05','07','09'};

num_epoch = 100;
dropout_rate = 0.0;


hydro_model = {'hymod','sacsma','vic'};


for hhh = 1
hydro_sel = hydro_model{hhh};
for kk=1:10
    
for www = 3


trial = kk;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       LOAD CDEC BASIN DATA (14 basins)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;

%----------------------------- Livneh CLIMATE ------------------------------------
basin_list = {'SHA','BND','ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB'};
basin_list_num = (1:length(basin_list));

% SPLIT DATA INTO TRAIN (1988-2003) & TEST (2004-2013)
syr_train = 1988;
eyr_train = 2003;
syr_test = 2004;
eyr_test = 2013;

% Livneh & PGF CLIMATE for each basin
datemat_livneh = datevec(datenum([1950 1 1]):datenum([2013 12 31]));
datemat_livneh_train = datevec(datenum([1988 1 1]):datenum([2003 12 31]));

basin_pr_train_all = nan(size(datemat_livneh_train,1),length(basin_list));
basin_tmax_train_all = nan(size(datemat_livneh_train,1),length(basin_list));
basin_tmin_train_all = nan(size(datemat_livneh_train,1),length(basin_list));
basin_daylen_train_all = nan(size(datemat_livneh_train,1),length(basin_list));
for i = 1:length(basin_list)
    
    basin_id = basin_list{i};
    
    basin_climate_livneh = load(['../meteo/livneh/fixed/meteo_livneh_CDEC_watershed_',basin_id,'.txt']);
    basin_daylen = load(['../meteo/LenDayLight/daylight_CDEC_watershed_',basin_id,'.txt']);
    % basin_climate : Year Mon Day PRCP(mm/day) Tmax(C) Tmin(C) Wind(m/s)
    basin_pr = basin_climate_livneh(:,4);
    basin_tmax = basin_climate_livneh(:,5);
    basin_tmin = basin_climate_livneh(:,6);
    basin_daylen = basin_daylen(:,4);

    
    sind_train = find(datemat_livneh(:,1)==syr_train & datemat_livneh(:,2)==1 & datemat_livneh(:,3)==1);
    eind_train = find(datemat_livneh(:,1)==eyr_train & datemat_livneh(:,2)==12 & datemat_livneh(:,3)==31);
    sind_test = find(datemat_livneh(:,1)==syr_test & datemat_livneh(:,2)==1 & datemat_livneh(:,3)==1);
    eind_test = find(datemat_livneh(:,1)==eyr_test & datemat_livneh(:,2)==12 & datemat_livneh(:,3)==31);         
    
    basin_pr_train = basin_pr(sind_train:eind_train);
    basin_tmax_train = basin_tmax(sind_train:eind_train);
    basin_tmin_train = basin_tmin(sind_train:eind_train);
    basin_daylen_train = basin_daylen(sind_train:eind_train);

    
    basin_pr_test = basin_pr(sind_test:eind_test);
    basin_tmax_test = basin_tmax(sind_test:eind_test);
    basin_tmin_test = basin_tmin(sind_test:eind_test);
    basin_daylen_test = basin_daylen(sind_test:eind_test);


    
    basin_pr_train_all(:,i) = basin_pr_train;
    basin_tmax_train_all(:,i) = basin_tmax_train;
    basin_tmin_train_all(:,i) = basin_tmin_train;
    basin_daylen_train_all(:,i) = basin_daylen_train;


    
    eval(['pr_train_',basin_id,' = basin_pr_train;'])
    eval(['tmax_train_',basin_id,' = basin_tmax_train;'])
    eval(['tmin_train_',basin_id,' = basin_tmin_train;'])
    eval(['daylen_train_',basin_id,' = basin_daylen_train;'])

    
    eval(['pr_test_',basin_id,' = basin_pr_test;'])
    eval(['tmax_test_',basin_id,' = basin_tmax_test;'])
    eval(['tmin_test_',basin_id,' = basin_tmin_test;'])  
    eval(['daylen_test_',basin_id,' = basin_daylen_test;'])
    
end

basin_pr_train_all_mean = mean(basin_pr_train_all(:));
basin_pr_train_all_std = std(basin_pr_train_all(:),1);
basin_tmax_train_all_mean = mean(basin_tmax_train_all(:));
basin_tmax_train_all_std = std(basin_tmax_train_all(:),1);
basin_tmin_train_all_mean = mean(basin_tmin_train_all(:));
basin_tmin_train_all_std = std(basin_tmin_train_all(:),1);
basin_daylen_train_all_mean = mean(basin_daylen_train_all(:));
basin_daylen_train_all_std = std(basin_daylen_train_all(:),1);

clear('basin_pr_train_all','basin_tmax_train_all','basin_tmin_train_all','basin_daylen_train_all')


% Standardizing Livneh CLIMATE for each basin
% e.g., pr_standard_train_SHA; pr_standard_test_SHA; tmax_standard_train_SHA; tmax_standard_test_SHA
for i = 1:length(basin_list)
    eval(['basin_pr_train = pr_train_',basin_list{i},';'])
    eval(['pr_standard_train_',basin_list{i},' = (basin_pr_train-basin_pr_train_all_mean)/basin_pr_train_all_std;'])
    %eval(['clear(''pr_train_',basin_list{i},''')'])
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
    
    eval(['basin_daylen_train = daylen_train_',basin_list{i},';'])
    eval(['daylen_standard_train_',basin_list{i},' = (basin_daylen_train-basin_daylen_train_all_mean)/basin_daylen_train_all_std;'])
    eval(['clear(''daylen_train_',basin_list{i},''')'])
    eval(['basin_daylen_test = daylen_test_',basin_list{i},';'])
    eval(['daylen_standard_test_',basin_list{i},' = (basin_daylen_test-basin_daylen_train_all_mean)/basin_daylen_train_all_std;'])
    eval(['clear(''daylen_test_',basin_list{i},''')'])
    
    
end

%------------------------ STATIC ATTRIBUTE --------------------------------
cdec_basin_static_feature = xlsread('../static/CatchmentStaticFeature.xlsx');
static_Area = cdec_basin_static_feature(1,:);
static_Lat = cdec_basin_static_feature(2,:);
static_Long = cdec_basin_static_feature(3,:);
static_Elev = cdec_basin_static_feature(4,:);
static_ElevRange = cdec_basin_static_feature(5,:);
static_Slp = cdec_basin_static_feature(6,:);
static_SlpRange = cdec_basin_static_feature(7,:);
static_Aspect = cdec_basin_static_feature(8,:);
static_FlowLen = cdec_basin_static_feature(9,:);
static_WaterFrac = cdec_basin_static_feature(10,:);
static_BareFrac = cdec_basin_static_feature(11,:);
static_UrbanFrac = cdec_basin_static_feature(12,:);
static_VegeFrac_Everg = cdec_basin_static_feature(13,:);
static_VegeFrac_Decidu = cdec_basin_static_feature(14,:);
static_VegeFrac_Mix = cdec_basin_static_feature(15,:);
static_VegeFrac_Wood = cdec_basin_static_feature(16,:);
static_VegeFrac_WoodGrass = cdec_basin_static_feature(17,:);
static_VegeFrac_CloseShrub = cdec_basin_static_feature(18,:);
static_VegeFrac_OpenShrub = cdec_basin_static_feature(19,:);
static_VegeFrac_Grass = cdec_basin_static_feature(20,:);
static_VegeFrac_Crop = cdec_basin_static_feature(21,:);
static_SoilDep_Pelletier = cdec_basin_static_feature(22,:);
static_SoilDep_Statsgo = cdec_basin_static_feature(23,:);
static_SoilFrac_Clay = cdec_basin_static_feature(24,:);
static_SoilFrac_ClayLoam = cdec_basin_static_feature(25,:);
static_SoilFrac_Loam = cdec_basin_static_feature(26,:);
static_SoilFrac_LoamSand = cdec_basin_static_feature(27,:);
static_SoilFrac_Sand = cdec_basin_static_feature(28,:);
static_SoilFrac_SandClayLoam = cdec_basin_static_feature(29,:);
static_SoilFrac_SiltClay = cdec_basin_static_feature(30,:);
static_SoilFrac_SiltClayLoam = cdec_basin_static_feature(31,:);
static_SoilFrac_SiltLoam = cdec_basin_static_feature(32,:);
static_SoilFrac_SandLoam = cdec_basin_static_feature(33,:);
static_SoilFrac_Other = cdec_basin_static_feature(34,:);

static_Area_stdr = (static_Area-mean(static_Area))/std(static_Area,1);
static_Lat_stdr = (static_Lat-mean(static_Lat))/std(static_Lat,1);
static_Long_stdr = (static_Long-mean(static_Long))/std(static_Long,1);
static_Elev_stdr = (static_Elev-mean(static_Elev))/std(static_Elev,1);
static_ElevRange_stdr = (static_ElevRange-mean(static_ElevRange))/std(static_ElevRange,1);
static_Slp_stdr = (static_Slp-mean(static_Slp))/std(static_Slp,1);
static_SlpRange_stdr = (static_SlpRange-mean(static_SlpRange))/std(static_SlpRange,1);
static_Aspect_stdr = (static_Aspect-mean(static_Aspect))/std(static_Aspect,1);
static_FlowLen_stdr = (static_FlowLen-mean(static_FlowLen))/std(static_FlowLen,1);
static_WaterFrac_stdr = (static_WaterFrac-mean(static_WaterFrac))/std(static_WaterFrac,1);
static_BareFrac_stdr = (static_BareFrac-mean(static_BareFrac))/std(static_BareFrac,1);
static_UrbanFrac_stdr = (static_UrbanFrac-mean(static_UrbanFrac))/std(static_UrbanFrac,1);
static_VegeFrac_Everg_stdr = (static_VegeFrac_Everg-mean(static_VegeFrac_Everg))/std(static_VegeFrac_Everg,1);
static_VegeFrac_Decidu_stdr = (static_VegeFrac_Decidu-mean(static_VegeFrac_Decidu))/std(static_VegeFrac_Decidu,1);
static_VegeFrac_Mix_stdr = (static_VegeFrac_Mix-mean(static_VegeFrac_Mix))/std(static_VegeFrac_Mix,1);
static_VegeFrac_Wood_stdr = (static_VegeFrac_Wood-mean(static_VegeFrac_Wood))/std(static_VegeFrac_Wood,1);
static_VegeFrac_WoodGrass_stdr = (static_VegeFrac_WoodGrass-mean(static_VegeFrac_WoodGrass))/std(static_VegeFrac_WoodGrass,1);
static_VegeFrac_CloseShrub_stdr = (static_VegeFrac_CloseShrub-mean(static_VegeFrac_CloseShrub))/std(static_VegeFrac_CloseShrub,1);
static_VegeFrac_OpenShrub_stdr = (static_VegeFrac_OpenShrub-mean(static_VegeFrac_OpenShrub))/std(static_VegeFrac_OpenShrub,1);
static_VegeFrac_Grass_stdr = (static_VegeFrac_Grass-mean(static_VegeFrac_Grass))/std(static_VegeFrac_Grass,1);
static_VegeFrac_Crop_stdr = (static_VegeFrac_Crop-mean(static_VegeFrac_Crop))/std(static_VegeFrac_Crop,1);
static_SoilDep_Pelletier_stdr = (static_SoilDep_Pelletier-mean(static_SoilDep_Pelletier))/std(static_SoilDep_Pelletier,1);
static_SoilDep_Statsgo_stdr = (static_SoilDep_Statsgo-mean(static_SoilDep_Statsgo))/std(static_SoilDep_Statsgo,1);
static_SoilFrac_Clay_stdr = (static_SoilFrac_Clay-mean(static_SoilFrac_Clay))/std(static_SoilFrac_Clay,1);
static_SoilFrac_ClayLoam_stdr = (static_SoilFrac_ClayLoam-mean(static_SoilFrac_ClayLoam))/std(static_SoilFrac_ClayLoam,1);
static_SoilFrac_Loam_stdr = (static_SoilFrac_Loam-mean(static_SoilFrac_Loam))/std(static_SoilFrac_Loam,1);
static_SoilFrac_LoamSand_stdr = (static_SoilFrac_LoamSand-mean(static_SoilFrac_LoamSand))/std(static_SoilFrac_LoamSand,1);
static_SoilFrac_Sand_stdr = (static_SoilFrac_Sand-mean(static_SoilFrac_Sand))/std(static_SoilFrac_Sand,1);
static_SoilFrac_SandClayLoam_stdr = (static_SoilFrac_SandClayLoam-mean(static_SoilFrac_SandClayLoam))/std(static_SoilFrac_SandClayLoam,1);
static_SoilFrac_SiltClay_stdr = (static_SoilFrac_SiltClay-mean(static_SoilFrac_SiltClay))/std(static_SoilFrac_SiltClay,1);
static_SoilFrac_SiltClayLoam_stdr = (static_SoilFrac_SiltClayLoam-mean(static_SoilFrac_SiltClayLoam))/std(static_SoilFrac_SiltClayLoam,1);
static_SoilFrac_SiltLoam_stdr = (static_SoilFrac_SiltLoam-mean(static_SoilFrac_SiltLoam))/std(static_SoilFrac_SiltLoam,1);
static_SoilFrac_SandLoam_stdr = (static_SoilFrac_SandLoam-mean(static_SoilFrac_SandLoam))/std(static_SoilFrac_SandLoam,1);
static_SoilFrac_Other_stdr = (static_SoilFrac_Other-mean(static_SoilFrac_Other))/std(static_SoilFrac_Other,1);


%------------------------ STREAMFLOW --------------------------------

% CDEC STREAMFLOW for each basin
datemat_train = datevec(datenum([1988 10 1]):datenum([2003 9 30]));
datemat_test = datevec(datenum([2003 10 1]):datenum([2013 9 30]));
basin_all_id_std_size = nan(length(basin_list),3);
for i = 1:length(basin_list)   
    
    basin_id = basin_list{i}; 
    basin_data = load(['../CDEC/FNF_mm_QC/FNF_',basin_id,'_mm.txt']); 
        
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
    basin_q_test(basin_q_test(:,4)<0,4) = nan;
    

    basin_all_id_std_size(i,:) = [i, nanstd(basin_q_train(:,4),1),sum(~isnan(basin_q_train(:,4)))];
    
    eval(['q_train_',basin_id,' = basin_q_train;'])
    eval(['q_test_',basin_id,' = basin_q_test;'])
    
end
%toc; % takes 20 sec from the very beginning


%------------------------ EVAPOTRANSPIRATION --------------------------------
% AET for each basin
datemat_train = datevec(datenum([1988 10 1]):datenum([2003 9 30]));
datemat_test = datevec(datenum([2003 10 1]):datenum([2013 9 30]));
%hydro_sel = 'hymod';
for i = 1:length(basin_list)   
    
    basin_id = basin_list{i}; 
    basin_data = load(['/home/ca-lstm/hydromodel/',hydro_sel,'/simaet_watershed/simaet_',hydro_sel,'_',basin_id,'.txt']); 
        
    sind_train = find(basin_data(:,1)==datemat_train(1,1) & basin_data(:,2)==datemat_train(1,2) & basin_data(:,3)==datemat_train(1,3));
    eind_train = find(basin_data(:,1)==datemat_train(end,1) & basin_data(:,2)==datemat_train(end,2) & basin_data(:,3)==datemat_train(end,3));
    sind_test = find(basin_data(:,1)==datemat_test(1,1) & basin_data(:,2)==datemat_test(1,2) & basin_data(:,3)==datemat_test(1,3));
    eind_test = find(basin_data(:,1)==datemat_test(end,1) & basin_data(:,2)==datemat_test(end,2) & basin_data(:,3)==datemat_test(end,3));  
    
    if isempty(sind_train)
        datenum_first = datenum(basin_data(1,1:3));
        datemat_miss = datevec(datenum(datemat_train(1,:)):datenum_first-1);
        basin_et_train = [[datemat_miss(:,1:3),repmat(-999,size(datemat_miss,1),1)];basin_data(1:eind_train,:)];
    else
        basin_et_train = basin_data(sind_train:eind_train,:);
    end
    if isempty(eind_test)
        datenum_last = datenum(basin_data(end,1:3));
        datemat_miss = datevec(datenum_last+1:datenum(datemat_test(end,:)));
        basin_et_test = [basin_data(sind_test:end,:);[datemat_miss(:,1:3),repmat(-999,size(datemat_miss,1),1)]];
    else
        basin_et_test = basin_data(sind_test:eind_test,:);
    end
    
    
    eval(['et_train_',basin_id,' = basin_et_train;'])
    eval(['et_test_',basin_id,' = basin_et_test;'])
    
end


%------------------------ DELTA STORAGE --------------------------------
% dS for each basin
datemat_train = datevec(datenum([1988 10 1]):datenum([2003 9 30]));
datemat_test = datevec(datenum([2003 10 1]):datenum([2013 9 30]));

for i = 1:length(basin_list)   
    
    basin_id = basin_list{i}; 
    basin_data = load(['/home/ca-lstm/hydromodel/',hydro_sel,'/simsm_watershed/simsm_',hydro_sel,'_',basin_id,'.txt']); 
        
    sind_train = find(basin_data(:,1)==datemat_train(1,1) & basin_data(:,2)==datemat_train(1,2) & basin_data(:,3)==datemat_train(1,3));
    eind_train = find(basin_data(:,1)==datemat_train(end,1) & basin_data(:,2)==datemat_train(end,2) & basin_data(:,3)==datemat_train(end,3));
    sind_test = find(basin_data(:,1)==datemat_test(1,1) & basin_data(:,2)==datemat_test(1,2) & basin_data(:,3)==datemat_test(1,3));
    eind_test = find(basin_data(:,1)==datemat_test(end,1) & basin_data(:,2)==datemat_test(end,2) & basin_data(:,3)==datemat_test(end,3));  
    
    if isempty(sind_train)
        datenum_first = datenum(basin_data(1,1:3));
        datemat_miss = datevec(datenum(datemat_train(1,:)):datenum_first-1);
        basin_ds_train = [[datemat_miss(:,1:3),repmat(-999,size(datemat_miss,1),1)];basin_data(1:eind_train,:)];
    else
        basin_ds_train = basin_data(sind_train:eind_train,:);
    end
    if isempty(eind_test)
        datenum_last = datenum(basin_data(end,1:3));
        datemat_miss = datevec(datenum_last+1:datenum(datemat_test(end,:)));
        basin_ds_test = [basin_data(sind_test:end,:);[datemat_miss(:,1:3),repmat(-999,size(datemat_miss,1),1)]];
    else
        basin_ds_test = basin_data(sind_test:eind_test,:);
    end
    
    
    eval(['ds_train_',basin_id,' = basin_ds_train;'])
    eval(['ds_test_',basin_id,' = basin_ds_test;'])
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    LSTM INPUT FOR TRAINING & TEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
num_feature_dynamic = 4;
seq_length = 270;
tot_sample_num = sum(basin_all_id_std_size(:,end));
XTrain_basin_id = nan(tot_sample_num,1);
XTrain = cell(tot_sample_num,1);
YTrain = cell(tot_sample_num,1);
YTrain_datenum = cell(tot_sample_num,1);
n=1;
for i = 1:length(basin_list)
    
    basin_id = basin_list{i};
    basin_id_num = basin_list_num(i);
    
    eval(['q_basin = q_train_',basin_id,';'])
    eval(['et_basin = et_train_',basin_id,';'])
    eval(['ds_basin = ds_train_',basin_id,';'])
    eval(['pr_basin_raw = pr_train_',basin_id,';'])
    
    eval(['pr_basin = pr_standard_train_',basin_id,';'])
    eval(['tmax_basin = tmax_standard_train_',basin_id,';'])
    eval(['tmin_basin = tmin_standard_train_',basin_id,';'])   
    eval(['daylen_basin = daylen_standard_train_',basin_id,';'])


    for j = 1:size(q_basin,1)
        
        if ~isnan(q_basin(j,end))
                        
            meteo_ind = find(datemat_livneh_train(:,1)==q_basin(j,1)&datemat_livneh_train(:,2)==q_basin(j,2)&datemat_livneh_train(:,3)==q_basin(j,3));
            x_seq = nan(num_feature_dynamic,seq_length);
            nn = 1;
            for tt = meteo_ind-seq_length+1:meteo_ind
                x_seq(:,nn) = [pr_basin(tt);tmax_basin(tt);tmin_basin(tt);daylen_basin(tt);...                   
                    ];
                nn=nn+1;
            end
            XTrain{n} = x_seq;
            YTrain{n} = [q_basin(j,4);et_basin(j,4);ds_basin(j,4)];
            YTrain_datenum{n} = datenum(q_basin(j,1:3));
            XTrain_basin_id(n) = basin_id_num;
            n=n+1;
        end       
    end
    
end
%toc; % takes about 10 min


for i = 1:length(basin_list)   
    basin_id = basin_list{i};
    eval(['clear(''pr_standard_train_',basin_id,''')'])
    eval(['clear(''tmax_standard_train_',basin_id,''')'])
    eval(['clear(''tmin_standard_train_',basin_id,''')'])    
    eval(['clear(''daylen_standard_train_',basin_id,''')'])

end
% toc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             LSTM TRAINING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- Loss Function Type -------------------------
loss_type = 'mse'; % either 'mse' or 'nse'

%--------------------- Hyperparameter being tuned -------------------------
%num_epoch = 50; % number of epoch

%------------------- Hyperparameter not being tuned -----------------------
% Number of input features
num_feature_static = 34;
num_input_lyr_unit = num_feature_dynamic+num_feature_static;

% Mini-batch size for SGD
batch_size = 1024;
% Adam optimizer learning rate 
alpha = .001; 
B1 = 0.9; 
B2 = 0.999; 
e = 10^-8; 
% Number of hidden layer units
num_hiddn_lyr_unit_1 = 128;
% Number of output layer units
num_output_lyr_unit = 3;
% Dropout rate
%dropout_rate = 0.10; % 0 is no dropout node; 0.2: on average 20% of hidden nodes are dropped out

%------------------- LSTM Weight Initialization ---------------------------
% Dimensions of Input, Hidden, Output Weights 
U_dim = [num_hiddn_lyr_unit_1,num_input_lyr_unit];
W_dim = [num_hiddn_lyr_unit_1,num_hiddn_lyr_unit_1];
V_dim = [num_output_lyr_unit,num_hiddn_lyr_unit_1];
% Row indies for each gate weights
iInd = (1:num_hiddn_lyr_unit_1*1);
fInd = (num_hiddn_lyr_unit_1*1+1:num_hiddn_lyr_unit_1*2);
oInd = (num_hiddn_lyr_unit_1*2+1:num_hiddn_lyr_unit_1*3);
gInd = (num_hiddn_lyr_unit_1*3+1:num_hiddn_lyr_unit_1*4);
ifoInd = [iInd fInd oInd];
% Xavier Glorot Initialization - Input Weight
r_U = sqrt(6/sum(U_dim));
U = -r_U+2*r_U*rand(U_dim(1)*4,U_dim(2),'single','gpuArray');
% Xavier Glorot Initialization - Hidden Weight
r_W = sqrt(6/sum(W_dim));
W = -r_W+2*r_W*rand(W_dim(1)*4,W_dim(2),'single','gpuArray');
% Xavier Glorot Initialization - Output Weight
r_V = sqrt(6/sum(V_dim));
V = -r_V+2*r_V*rand(V_dim,'single','gpuArray');
% Bias Initialization
b = ones(U_dim(1)*4,1,'single','gpuArray');
bv = ones(V_dim(1),1,'single','gpuArray');


%---------------- Momentum & AdaGrad Initialization for Adam --------------
V_mom_0 = zeros(V_dim,'single','gpuArray');
bv_mom_0 = zeros(V_dim(1),1,'single','gpuArray');
U_mom_0 = zeros(U_dim(1)*4,U_dim(2),'single','gpuArray');
W_mom_0 = zeros(W_dim(1)*4,W_dim(2),'single','gpuArray');
b_mom_0 = zeros(U_dim(1)*4,1,'single');

V_adagrad_0 = zeros(V_dim,'single','gpuArray');
bv_adagrad_0 = zeros(V_dim(1),1,'single','gpuArray');
W_adagrad_0 = zeros(W_dim(1)*4,W_dim(2),'single','gpuArray');
U_adagrad_0 = zeros(U_dim(1)*4,U_dim(2),'single','gpuArray');
b_adagrad_0 = zeros(U_dim(1)*4,1,'single','gpuArray');

%--------------------------- LSTM Training -------------------------------- 
XTrain_calib = XTrain;
YTrain_calib = YTrain;
YTrain_datenum_calib = YTrain_datenum;
XTrain_calib_basin_id = XTrain_basin_id;

% Number of iterations based on mini batch size
num_iter = max([1,round(length(XTrain_calib)/batch_size)]);
batch_size_vec = repmat(min([batch_size,length(XTrain_calib)]),num_iter,1);
extra = length(XTrain_calib) - min([batch_size,length(XTrain_calib)])*num_iter;
if extra > 0
    n = 1;
    while extra > 0 
        batch_size_vec(n) = batch_size_vec(n)+1;
        extra = extra - 1;
        n = n+1;
        if n > num_iter; n=1; end
    end    
elseif extra < 0 
    n = 1;
    while extra < 0 
        batch_size_vec(n) = batch_size_vec(n)-1;
        extra = extra + 1;
        n = n+1;
        if n > num_iter; n=1; end
    end
end
batch_ind = cumsum(batch_size_vec)-batch_size_vec+1;


% LSTM Forward & Back propagation 
batch_info_selbasin = zeros(max(batch_size_vec),num_epoch*num_iter,'single');
batch_info_datenum = zeros(max(batch_size_vec),num_epoch*num_iter,'single');
batch_info_Yhat = zeros(max(batch_size_vec),num_epoch*num_iter,'single');
batch_info_Yobs = zeros(max(batch_size_vec),num_epoch*num_iter,'single');

tic;
loss_calib = nan(num_epoch*num_iter,1);
loss_calib_1st = nan(num_epoch*num_iter,1);
loss_calib_2nd = nan(num_epoch*num_iter,1);
loss_calib_3rd = nan(num_epoch*num_iter,1);
n = 1;
for i = 1:num_epoch
    
    % Shuffle Mini Batch Every Epoch
    p_batch = randperm(size(XTrain_calib,1));
    XTrain_calib_shuffle = XTrain_calib(p_batch); 
    YTrain_calib_shuffle = YTrain_calib(p_batch);
    YTrain_datenum_calib_shuffle = YTrain_datenum_calib(p_batch,:);
    XTrain_calib_basin_id_shuffle = XTrain_calib_basin_id(p_batch);
    
    for j = 1:num_iter
% tic;      
        XTrain_mini = XTrain_calib_shuffle(batch_ind(j):batch_ind(j)+batch_size_vec(j)-1);
        YTrain_mini = YTrain_calib_shuffle(batch_ind(j):batch_ind(j)+batch_size_vec(j)-1);
        YTrain_datenum_mini = YTrain_datenum_calib_shuffle(batch_ind(j):batch_ind(j)+batch_size_vec(j)-1);
        XTrain_calib_basin_id_mini = XTrain_calib_basin_id_shuffle(batch_ind(j):batch_ind(j)+batch_size_vec(j)-1);
        minibatch_size = size(XTrain_mini,1);
        
        if strcmp(loss_type,'nse')           
            Y_std = basin_all_id_std_size(XTrain_calib_basin_id_mini,2);
            err_scale = single(1./(basin_all_id_std_size(XTrain_calib_basin_id_mini,3).*(Y_std.^2)*size(basin_all_id_std_size,1)));
        elseif strcmp(loss_type,'mse')
            err_scale = ones(size(XTrain_mini),'single');
        end
        
        attr_static = zeros(num_feature_static,minibatch_size);
        
        attr_static(1,:) = static_Area_stdr(XTrain_calib_basin_id_mini);
        attr_static(2,:) = static_Lat_stdr(XTrain_calib_basin_id_mini);
        attr_static(3,:) = static_Long_stdr(XTrain_calib_basin_id_mini);
        attr_static(4,:) = static_Elev_stdr(XTrain_calib_basin_id_mini);
        attr_static(5,:) = static_ElevRange_stdr(XTrain_calib_basin_id_mini);
        attr_static(6,:) = static_Slp_stdr(XTrain_calib_basin_id_mini);
        attr_static(7,:) = static_SlpRange_stdr(XTrain_calib_basin_id_mini);
        attr_static(8,:) = static_Aspect_stdr(XTrain_calib_basin_id_mini);
        attr_static(9,:) = static_FlowLen_stdr(XTrain_calib_basin_id_mini);
        attr_static(10,:) = static_WaterFrac_stdr(XTrain_calib_basin_id_mini);
        attr_static(11,:) = static_BareFrac_stdr(XTrain_calib_basin_id_mini);
        attr_static(12,:) = static_UrbanFrac_stdr(XTrain_calib_basin_id_mini);
        attr_static(13,:) = static_VegeFrac_Everg_stdr(XTrain_calib_basin_id_mini);
        attr_static(14,:) = static_VegeFrac_Decidu_stdr(XTrain_calib_basin_id_mini);
        attr_static(15,:) = static_VegeFrac_Mix_stdr(XTrain_calib_basin_id_mini);
        attr_static(16,:) = static_VegeFrac_Wood_stdr(XTrain_calib_basin_id_mini);
        attr_static(17,:) = static_VegeFrac_WoodGrass_stdr(XTrain_calib_basin_id_mini);
        attr_static(18,:) = static_VegeFrac_CloseShrub_stdr(XTrain_calib_basin_id_mini);
        attr_static(19,:) = static_VegeFrac_OpenShrub_stdr(XTrain_calib_basin_id_mini);
        attr_static(20,:) = static_VegeFrac_Grass_stdr(XTrain_calib_basin_id_mini);
        attr_static(21,:) = static_VegeFrac_Crop_stdr(XTrain_calib_basin_id_mini);
        attr_static(22,:) = static_SoilDep_Pelletier_stdr(XTrain_calib_basin_id_mini);
        attr_static(23,:) = static_SoilDep_Statsgo_stdr(XTrain_calib_basin_id_mini);
        attr_static(24,:) = static_SoilFrac_Clay_stdr(XTrain_calib_basin_id_mini);
        attr_static(25,:) = static_SoilFrac_ClayLoam_stdr(XTrain_calib_basin_id_mini);
        attr_static(26,:) = static_SoilFrac_Loam_stdr(XTrain_calib_basin_id_mini);
        attr_static(27,:) = static_SoilFrac_LoamSand_stdr(XTrain_calib_basin_id_mini);
        attr_static(28,:) = static_SoilFrac_Sand_stdr(XTrain_calib_basin_id_mini);
        attr_static(29,:) = static_SoilFrac_SandClayLoam_stdr(XTrain_calib_basin_id_mini);
        attr_static(30,:) = static_SoilFrac_SiltClay_stdr(XTrain_calib_basin_id_mini);
        attr_static(31,:) = static_SoilFrac_SiltClayLoam_stdr(XTrain_calib_basin_id_mini);
        attr_static(32,:) = static_SoilFrac_SiltLoam_stdr(XTrain_calib_basin_id_mini);
        attr_static(33,:) = static_SoilFrac_SandLoam_stdr(XTrain_calib_basin_id_mini);
        attr_static(34,:) = static_SoilFrac_Other_stdr(XTrain_calib_basin_id_mini);

        X_static = reshape(repmat(attr_static(:),1,seq_length),num_feature_static,minibatch_size,seq_length);
        X_dynamic = permute(reshape([XTrain_mini{:}],num_feature_dynamic,seq_length,minibatch_size),[1,3,2]);
        
        X_all = zeros(num_input_lyr_unit,minibatch_size,seq_length,'single');
        X_all(1:num_feature_dynamic,:,:) = X_dynamic;
        X_all(num_feature_dynamic+1:num_input_lyr_unit,:,:) = X_static; 
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LSTM FORWARD %%%%%%%%%%%%%%%%%%%%%%%
%         X = gpuArray(X_all);
%         H_all = zeros(num_hiddn_lyr_unit_1,minibatch_size,seq_length,'single');
%         C_all = zeros(num_hiddn_lyr_unit_1,minibatch_size,seq_length,'single');      
%         C = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');  
%         H = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');        
%         
%         for tt = 1:seq_length 
%             % Gate operation    
%             G = U*X(:,:,tt) + W*H + b;
%             G(ifoInd,:) = 1./(1+exp(-G(ifoInd,:)));
%             G(gInd,:)=tanh(G(gInd,:));
%                
%             % Cell state update
%             C = G(gInd,:).*G(iInd,:)+G(fInd,:).*C;
%             C_all(:,:,tt) = C;
% 
%             % LSTM layer output
%             H = tanh(C).*G(oInd, :);
%             H_all(:,:,tt) = H;
% 
%         end                
%         
%         % DROPOUT for LSTM hidden node output
%         dropInd = rand(num_hiddn_lyr_unit_1,1) < dropout_rate;
%         dropoutScaleFactor = 1-dropout_rate;
%         H = H/dropoutScaleFactor;
%         H(dropInd,:) = 0;
%        
%         % OUTPUT LAYER
%         Yhat = V*H + bv;
%         Yhat(1,Yhat(1,:)<0) = 0; % ReLU activation function
%         Yobs = single(reshape([YTrain_mini{:}],3,minibatch_size));

        X = gpuArray(X_all);
        Yobs = single(gpuArray(reshape([YTrain_mini{:}],3,minibatch_size)));
        h0 = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');
        c0 = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');  
        dropInd = zeros(num_hiddn_lyr_unit_1,1,'single','gpuArray');
        r = rand(num_hiddn_lyr_unit_1,1);
        for k = 1:num_hiddn_lyr_unit_1
            if r(k) < dropout_rate
                dropInd(k) = 1;
            end
        end
      
        % Executing CUDA Kernel function for LSTM Forward Process
        [Yhat, Loss, H_all, C_all] = lstm_forward_mexcuda(...
            U(iInd,:),U(fInd,:),U(oInd,:),U(gInd,:),...
            W(iInd,:),W(fInd,:),W(oInd,:),W(gInd,:),...
            b(iInd,:),b(fInd,:),b(oInd,:),b(gInd,:),...
            V,bv,...
            X,h0,c0,Yobs,gpuArray(err_scale),...
            single(gpuArray(dropout_rate)), dropInd...
            );


        Yhat(1,(Yhat(1,:)<0)) = 0; % ReLU activation function
        
        
        L_mse = 1/2*(Yhat(1,:)-Yobs(1,:)).^2;        
        L_mse_batch = sum(mean(L_mse,1))/minibatch_size;  
               
        L_mse_et = 1/2*(Yhat(2,:)-Yobs(2,:)).^2;        
        L_mse_et_batch = sum(mean(L_mse_et,1))/minibatch_size;
        
        L_mse_ds = 1/2*(Yhat(3,:)-Yobs(3,:)).^2;        
        L_mse_ds_batch = sum(mean(L_mse_ds,1))/minibatch_size;
        


        %Loss = w_q(www)*L_mse_batch+w_et(www)*L_mse_et_batch;
        Loss = L_mse_batch + L_mse_et_batch + L_mse_ds_batch;
        Loss_1st = L_mse_batch;
        Loss_2nd = L_mse_et_batch;
        Loss_3rd = L_mse_ds_batch;
             

        
  
% % toc;        
% % tic;       
%         % FORWARD
%         X = gpuArray(X_all);
%         Yobs = single(gpuArray([YTrain_mini{:}]));
%         h0 = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');
%         c0 = zeros(num_hiddn_lyr_unit_1,minibatch_size,'single','gpuArray');  
%         dropInd = zeros(num_hiddn_lyr_unit_1,1,'single','gpuArray');
%         r = rand(num_hiddn_lyr_unit_1,1);
%         for k = 1:num_hiddn_lyr_unit_1
%             if r(k) < dropout_rate
%                 dropInd(k) = 1;
%             end
%         end
% % tic;      
%         % Executing CUDA Kernel function for LSTM Forward Process
%         [Yhat, Loss, H_all, C_all] = lstm_forward_mexcuda(...
%             U(iInd,:),U(fInd,:),U(oInd,:),U(gInd,:),...
%             W(iInd,:),W(fInd,:),W(oInd,:),W(gInd,:),...
%             b(iInd,:),b(fInd,:),b(oInd,:),b(gInd,:),...
%             V,bv,...
%             X,h0,c0,Yobs,gpuArray(err_scale),...
%             single(gpuArray(dropout_rate)), dropInd...
%             );
% % wait(gpu); toc;   
% % tic;
        
        batch_info_selbasin(1:minibatch_size,n) = XTrain_calib_basin_id_mini;
        batch_info_datenum(1:minibatch_size,n) = [YTrain_datenum_mini{:}];
        batch_info_Yhat(1:minibatch_size,n) = Yhat(1,:);
        batch_info_Yobs(1:minibatch_size,n) = Yobs(1,:);

        % BACKWARD  
        H_end = H_all(:,:,end)/(1-dropout_rate);
        H_end(dropInd==1,:) = 0;
        
        d_relu = zeros(size(Yhat(1,:)));
        d_relu(Yhat(1,:)>0) = 1;
        
 
        %dy_q = ((Yhat(1,:)-Yobs(1,:)).*d_relu+(Yhat(1,:)+Yhat(2,:)-(Yobs(2,:)+Yobs(1,:))).*d_relu)* w_mse(www); % WRONG
        dy_q = (Yhat(1,:)-Yobs(1,:)).*d_relu;
        dy_et = (Yhat(2,:)-Yobs(2,:));
        dy_ds = (Yhat(3,:)-Yobs(3,:));
        
        dV = [dy_q;dy_et;dy_ds]*H_end';
        dbv = [sum(dy_q);sum(dy_et);sum(dy_ds)];


        %%%%%%%%% Backpropatation through time
        
        %dX = zeros(size(X,1),minibatch_size,'gpuArray');
        dG = zeros(W_dim(1)*4,minibatch_size,'single','gpuArray');      
        dU = zeros(U_dim(1)*4,U_dim(2),'single','gpuArray');
        dW = zeros(W_dim(1)*4,W_dim(2),'single','gpuArray');
        db = zeros(U_dim(1)*4,1,'single','gpuArray');
        
        tt = seq_length;
        
        % Layer output derivative
        dH = V'*(Yhat-Yobs);
        
        % Tanh activation of the cell state
        tanhC = tanh(C_all(:,:,tt));
        
        % Determine the gates
        Gi = 1./(1+exp(-(U(iInd,:)*X(:,:,tt) + W(iInd,:)*H_all(:,:,tt-1) + b(iInd))));
        Gf = 1./(1+exp(-(U(fInd,:)*X(:,:,tt) + W(fInd,:)*H_all(:,:,tt-1) + b(fInd))));
        Go = 1./(1+exp(-(U(oInd,:)*X(:,:,tt) + W(oInd,:)*H_all(:,:,tt-1) + b(oInd))));
        Gg = tanh(U(gInd,:)*X(:,:,tt) + W(gInd,:)*H_all(:,:,tt-1) + b(gInd));

        
        % Cell state
        %C = (C - Gg.*Gi)./Gf;                
        
        % Output gate derivative
        dG(oInd,:) = (dH.*tanhC).*Go.*(1-Go);
        
        % Cell state derivative
        dC = dH.*Go.*(1-tanhC.^2);
        
        % Forget gate derivative
        %dG(fInd,:) = (dC.*C).*G(fInd,:).*(1-G(fInd,:));
        dG(fInd,:) = (dC.*C_all(:,:,tt-1)).*Gf.*(1-Gf);
        
        % Input gate derivative
        dG(iInd,:) = (dC.*Gg).*Gi.*(1-Gi);
        
        % Candidate derivative
        dG(gInd,:) = (dC.*Gi).*(1-(Gg).^2);
        
        % Input data derivative
        % dX(:, :, tt) = U(iInd,:)'*dG(iInd,:) + U(fInd,:)'*dG(fInd,:) + U(oInd,:)'*dG(oInd,:) + U(gInd,:)'*dG(gInd,:);
        
        % Input weights (U) derivative
        dU = dU + dG*X(:,:,tt)';        

        % Recurrent weights (W) derivative
        dW = dW + dG*H_all(:,:,tt-1)';
        
        % Bias derivative
        db = db + sum(dG, 2);
        
        for tt = seq_length-1:-1:2
            % Layer output derivative
            dH = W(iInd,:)'*dG(iInd,:) + W(fInd,:)'*dG(fInd,:) + W(oInd,:)'*dG(oInd,:) + W(gInd,:)'*dG(gInd,:);
            
            % Tanh activation of the cell state
            tanhC = tanh(C_all(:,:,tt));
            
            % Store Gf for tt+1
            Gfttp = Gf;
            
            % Determine the gates
            Gi = 1./(1+exp(-(U(iInd,:)*X(:,:,tt) + W(iInd,:)*H_all(:,:,tt-1) + b(iInd))));
            Gf = 1./(1+exp(-(U(fInd,:)*X(:,:,tt) + W(fInd,:)*H_all(:,:,tt-1) + b(fInd))));
            Go = 1./(1+exp(-(U(oInd,:)*X(:,:,tt) + W(oInd,:)*H_all(:,:,tt-1) + b(oInd))));
            Gg = tanh(U(gInd,:)*X(:,:,tt) + W(gInd,:)*H_all(:,:,tt-1) + b(gInd));
            
            
            % Cell state
            %C = (C - Gg.*Gi)./Gf; 
            
            % Output gate derivative
            dG(oInd,:) = (dH.*tanhC).*Go.*(1-Go);
            
            % Cell state derivative
            dC = dH.*Go.*(1-tanhC.^2) + dC.*Gfttp;           

            % Forget gate derivative
            dG(fInd,:) = (dC.*C_all(:,:,tt-1)).*Gf.*(1-Gf);

            % Input gate derivative
            dG(iInd,:) = (dC.*Gg).*Gi.*(1-Gi);

            % Candidate derivative
            dG(gInd,:) = (dC.*Gi).*(1-(Gg).^2);

            % Input data derivative
            % dX(:, :, tt) = U(iInd,:)'*dG(iInd,:) + U(fInd,:)'*dG(fInd,:) + U(oInd,:)'*dG(oInd,:) + U(gInd,:)'*dG(gInd,:);

            % Input weights (U) derivative
            dU = dU + dG*X(:,:,tt)';        

            % Recurrent weights (W) derivative
            dW = dW + dG*H_all(:,:,tt-1)';

            % Bias derivative
            db = db + sum(dG, 2);

        end   
        
        tt = 1;
        h0 = repmat(zeros(num_hiddn_lyr_unit_1,1,'single','gpuArray'),1,minibatch_size);
        c0 = zeros(num_hiddn_lyr_unit_1,1,'single','gpuArray');
        % Layer output derivative
        dH = W(iInd,:)'*dG(iInd,:) + W(fInd,:)'*dG(fInd,:) + W(oInd,:)'*dG(oInd,:) + W(gInd,:)'*dG(gInd,:);
        
        % Tanh activation of the cell state
        tanhC = tanh(C_all(:,:,tt));

        % Store Gf for tt+1
        Gfttp = Gf;
        
        % Determine the gates
        Gi = 1./(1+exp(-(U(iInd,:)*X(:,:,tt) + W(iInd,:)*h0 + b(iInd))));
        Gf = 1./(1+exp(-(U(fInd,:)*X(:,:,tt) + W(fInd,:)*h0 + b(fInd))));
        Go = 1./(1+exp(-(U(oInd,:)*X(:,:,tt) + W(oInd,:)*h0 + b(oInd))));
        Gg = tanh(U(gInd,:)*X(:,:,tt) + W(gInd,:)*h0 + b(gInd));
        
        % Output gate derivative
        dG(oInd,:) = (dH.*tanhC).*Go.*(1-Go);

        % Cell state derivative
        dC = dH.*Go.*(1-tanhC.^2) + dC.*Gfttp;           

        % Forget gate derivative
        dG(fInd,:) = (dC.*c0).*Gf.*(1-Gf);

        % Input gate derivative
        dG(iInd,:) = (dC.*Gg).*Gi.*(1-Gi);

        % Candidate derivative
        dG(gInd,:) = (dC.*Gi).*(1-(Gg).^2);

        % Input data derivative
        % dX(:, :, tt) = U(iInd,:)'*dG(iInd,:) + U(fInd,:)'*dG(fInd,:) + U(oInd,:)'*dG(oInd,:) + U(gInd,:)'*dG(gInd,:);

        % Input weights (U) derivative
        dU = dU + dG*X(:,:,tt)';        

        % Recurrent weights (W) derivative
        dW = dW + dG*h0';

        % Bias derivative
        db = db + sum(dG, 2);
% toc;
% tic;        
        
        V_mom  = B1*V_mom_0+(1-B1)*dV;
        bv_mom = B1*bv_mom_0+(1-B1)*dbv;
        W_mom = B1*W_mom_0+(1-B1)*dW;
        U_mom = B1*U_mom_0+(1-B1)*dU;
        b_mom = B1*b_mom_0+(1-B1)*db;
        
        V_mom_bc  = V_mom/(1-B1^n);
        bv_mom_bc = bv_mom/(1-B1^n);
        W_mom_bc = W_mom/(1-B1^n);
        U_mom_bc = U_mom/(1-B1^n);
        b_mom_bc = b_mom/(1-B1^n);
        
        V_adagrad = B2*V_adagrad_0+(1-B2)*(dV).^2;
        bv_adagrad = B2*bv_adagrad_0+(1-B2)*(dbv).^2;
        W_adagrad = B2*W_adagrad_0+(1-B2)*dW.^2;
        U_adagrad = B2*U_adagrad_0+(1-B2)*dU.^2;
        b_adagrad = B2*b_adagrad_0+(1-B2)*db.^2;
        
        V_adagrad_bc = V_adagrad/(1-B2^n);
        bv_adagrad_bc = bv_adagrad/(1-B2^n);
        W_adagrad_bc = W_adagrad/(1-B2^n);
        U_adagrad_bc = U_adagrad/(1-B2^n);
        b_adagrad_bc = b_adagrad/(1-B2^n);
        
        
        V = V-alpha*V_mom_bc./(sqrt(V_adagrad_bc)+e);
        bv = bv-alpha*bv_mom_bc./(sqrt(bv_adagrad_bc)+e);
        W = W-alpha*W_mom_bc./(sqrt(W_adagrad_bc)+e);
        U = U-alpha*U_mom_bc./(sqrt(U_adagrad_bc)+e);
        b = b-alpha*b_mom_bc./(sqrt(b_adagrad_bc)+e);        
        
        
        V_adagrad_0 = V_adagrad;
        bv_adagrad_0 = bv_adagrad;
        W_adagrad_0 = W_adagrad;
        U_adagrad_0 = U_adagrad;
        b_adagrad_0 = b_adagrad;       
        
        
        V_mom_0 = V_mom;
        bv_mom_0 = bv_mom;
        W_mom_0 = W_mom;
        U_mom_0 = U_mom;
        b_mom_0 = b_mom;
% toc;        
        
        if strcmp(loss_type,'nse')
            loss_calib(n) = 2*sum(mean(Loss,1));
            NSE = 1-loss_calib(n)*num_iter; % approx NSE
            elapse_t = toc;
            disp(['Epoch: ',num2str(i),'/',num2str(num_epoch),...
                '    Iteration: ',num2str(j),'/',num2str(num_iter),...
                '   NSE : ',num2str(NSE),...
                '   Loss : ',num2str(loss_calib(n)),...
                '   Elapsed Time: ',num2str(elapse_t), 'sec'])
        elseif strcmp(loss_type,'mse')
            %loss_calib(n) = sum(mean(Loss,1))/minibatch_size;
            %RMSE = sqrt(sum(mean(2*Loss,1))/minibatch_size);
            loss_calib(n) = Loss;
            loss_calib_1st(n) = Loss_1st;
            loss_calib_2nd(n) = Loss_2nd;
            loss_calib_3rd(n) = Loss_3rd;
            RMSE = sqrt(2*Loss);
            elapse_t = toc;
            disp(['Epoch: ',num2str(i),'/',num2str(num_epoch),...
                '    Iteration: ',num2str(j),'/',num2str(num_iter),...
                '   RMSE : ',num2str(RMSE),...
                '   Loss1 : ',num2str(loss_calib_1st(n)),...
                '   Loss2 : ',num2str(loss_calib_2nd(n)),...
                '   Loss3 : ',num2str(loss_calib_3rd(n)),...
                '   Elapsed Time: ',num2str(elapse_t), 'sec'])
        end        
        
        n = n+1;
        
    end
        
    
    
end


% toc;
% figure
% plot(loss_calib,'k'); hold on;
% plot(loss_calib_1st,'b'); hold on;
% plot(loss_calib_2nd,'r'); hold on;
% plot(loss_calib_3rd,'g'); hold on;
% xlabel('Number of Iteration')
% ylabel('Loss Function')
% title('CA-LSTM Loss Function Convergence')
% grid on;


grp_arr = nan(num_iter*num_epoch,1);
n=1;
for i=1:num_epoch
    grp_arr(n:n+num_iter-1) = repmat(i,num_iter,1);
    n=n+num_iter;
end

loss_calib_epoch = grpstats(loss_calib,grp_arr,'mean');
loss_calib_1st_epoch = grpstats(loss_calib_1st,grp_arr,'mean');
loss_calib_2nd_epoch = grpstats(loss_calib_2nd,grp_arr,'mean');
loss_calib_3rd_epoch = grpstats(loss_calib_3rd,grp_arr,'mean');

% figure
% plot(loss_calib_epoch,'k'); hold on;
% plot(loss_calib_1st_epoch,'b'); hold on;
% plot(loss_calib_2nd_epoch,'r'); hold on;
% plot(loss_calib_3rd_epoch,'g'); hold on;
% xlabel('Number of Epoch')
% ylabel('Loss Function')
% title('CA-LSTM Loss Function Convergence')
% grid on;

U_final = gather(U);
W_final = gather(W);
V_final = gather(V);
b_final = gather(b);
bv_final = gather(bv);

save(['hybrid_lstm_',hydro_sel,'_ETout_DSout_epoch',num2str(num_epoch),'_drop',num2str(dropout_rate*100),'_trial',num2str(trial),'.mat'],...
    'U_final','W_final','V_final','b_final','bv_final',...
    'loss_calib_epoch','loss_calib','loss_calib_1st_epoch','loss_calib_1st','loss_calib_2nd_epoch','loss_calib_2nd','loss_calib_3rd_epoch','loss_calib_3rd',...
    'batch_info_selbasin','batch_info_datenum','batch_info_Yhat','batch_info_Yobs');

end

end

end
