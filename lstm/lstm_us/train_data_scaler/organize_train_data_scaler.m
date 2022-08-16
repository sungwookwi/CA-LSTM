
clear
clc



for t = 1:10
    d = xlsread(['train_data_scaler_raw_t',num2str(t),'.csv']);
    
    dyn_mean = d(2,1:5);
    dyn_std = d(3,1:5);
    
    target_mean = d(4,1);
    target_std = d(5,1);
    
    attr_mean = d(6,1:26);
    attr_std = d(7,1:26);
    
    dyn_name = {'PRCP';'SRAD';'Tmax';'Tmin';'Vp'};
    dyn_save = [dyn_name,num2cell([dyn_mean',dyn_std'])];
    
    target_name = {'FLOW'};
    target_save = [target_name,num2cell([target_mean',target_std'])];
    

    attr_name = {'carbonate_rocks_frac';'geol_permeability';'frac_forest';'lai_max';'lai_diff';...
        'gvf_max';'gvf_diff';'p_mean';'pet_mean';'frac_snow';'aridity';'high_prec_freq';'high_prec_dur';...
        'low_prec_freq';'low_prec_dur';'elev_mean';'slope_mean';'area_gages2';'soil_depth_pelletier';...
        'soil_depth_statsgo';'soil_porosity';'soil_conductivity';'max_water_content';'sand_frac';'silt_frac';'clay_frac'};
    attr_save = [attr_name,num2cell([attr_mean',attr_std'])];
    
    xlswrite(['train_data_scaler_t',num2str(t),'.xls'],dyn_save,'dyn')
    xlswrite(['train_data_scaler_t',num2str(t),'.xls'],target_save,'target')
    xlswrite(['train_data_scaler_t',num2str(t),'.xls'],attr_save,'camel_attr')
    
end