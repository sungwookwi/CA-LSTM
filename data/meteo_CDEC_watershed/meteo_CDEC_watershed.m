
clear
clc


cdec_id = {'SHA','BND','ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB'};
tic;
for i = 2:length(cdec_id)
    hruinfo = load(['E:\StochasticWatershedModeling\SACSMA\hruinfo\HRUinfo_',cdec_id{i},'.txt']);
    
    datemat = datevec(datenum([1950 1 1]):datenum([2013 12 31]));
    hru_pr_tot = nan(size(datemat,1),size(hruinfo,1));
    hru_tas_tot = nan(size(datemat,1),size(hruinfo,1));
    for k = 1:size(hruinfo,1)
        hru_meteo = load(['E:\StochasticWatershedModeling\VIC\hru_meteo\hymod\meteo_',num2str(hruinfo(k,1),'%1.6f'),'_',num2str(hruinfo(k,2),'%1.6f')]);
        hru_pr_tot(:,k) = hru_meteo(:,4)*hruinfo(k,3)/sum(hruinfo(:,3));
        hru_tas_tot(:,k) = hru_meteo(:,5)*hruinfo(k,3)/sum(hruinfo(:,3));
    end
    
    hru_pr_basin = sum(hru_pr_tot,2);
    hru_tas_basin = sum(hru_tas_tot,2);
    
    
    savedata = [datemat(:,1:3),hru_pr_basin,hru_tas_basin];
    fid = fopen(['E:\StochasticWatershedModeling\meteo\meteo_CDEC_watershed\meteo_CDEC_watershed_',cdec_id{i},'.txt'],'w');
    fprintf(fid,'%4d %2d %2d %8.2f %8.2f\n',savedata');
    fclose(fid);
    toc;
end


