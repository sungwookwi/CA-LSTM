function [runoff, baseflow] = rout_lohmann(inflow_direct, inflow_base, flowlen, route_par, isOutlet)
%#codegen
% Routing model for Land Surface model outputs based on Lohmann routing
% model
% 
% Catchment's UH is represented by the Gamma distribution 
% River routing is based on the linearized Saint-Venant Eq 
% 
% Inputs
%	inflow_direct  =  Direct runoff from catchment (surface runoff, prompt subsurface runoff)
%	inflow_base    =  Baseflow from catchment (groundwater runoff, delayed subsurface runoff)
%   flowlen        =  Travel distance of water from catchment outlet to watershed outlet
%	route_par      =  UH & Lohmann routing parameters
%	isOutlet       =  Indicator telling whether the catchment is for watershed outlet 
%	
% Parameters
%	route_par(1)  =  Catchment's UH shape parameter (N)
%	route_par(2)  =  Catchment's UH scale parameter (K)
%	route_par(3)  =  Wave velocity in the linearized Saint-Venant equation(m/s) (VELO)
%	route_par(4)  =  Diffusivity in the linearized Saint-Venant equation(m2/s) (DIFF)
% 
% Outputs
%   runoff   = Runoff at the watershed outlet 
%   baseflow = Baseflow at the watershed outlet 
% 
% codegen codes
%	cfg = coder.config('mex');
%   cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%   codegen rout_lohmann -args {coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),0,zeros(1,4),0} -config cfg -report
% 
% Release Notes
% 	Written by Sungwook Wi(sungwookwi@gmail.com) 
%	Last Updated 05/11/2016 
% 

%------------------------------ Lohamann parameters -----------------------
N      = route_par(1);
K      = route_par(2);
VELO   = route_par(3);
DIFF   = route_par(4);
%--------------------------------------------------------------------------


%---- Base Time for HRU(watershed subunit) UH and channel routing UH ------
KE      =  12;          % Base time for HRU UH (day)
UH_DAY  =  96;          % Base time for river routing UH (day)
DT      =  3600;        % Time step in second for solving Saint-Venant equation. This will affect TMAX
TMAX    =  UH_DAY * 24; % Base time of river routing UH in hour because DT is for an hour
LE      =  48*50;       % Base time (hr) for Green function values
%--------------------------------------------------------------------------


%----- Derive Daily River Impulse Response Function(Green's function) -----
UH_river = zeros(1,UH_DAY);
if isOutlet == 1 % if the HUR contains the watershed outlet
    
    UH_river(1) =  1;
    
else % if not, calculate Green's function to solve the linearized Saint-Venant Eq
    
    t = 0;
    uhm_grid = zeros(LE,1);
    for k = 1:LE
        t = t + DT;
        
        pot = ((VELO*t-flowlen)^2)/(4*DIFF*t);
        if pot <= 69
            H = flowlen./(2*t.*sqrt(pi*t*DIFF)).*exp(-pot);
        else
            H = 0.0;
        end
        uhm_grid(k) = H;
        
    end
    
    if sum(uhm_grid) == 0
        uhm_grid(1) =  1.0;
    else
        uhm_grid = uhm_grid/sum(uhm_grid);
    end
        
    UHM      = uhm_grid;
    
    FR  = zeros(TMAX,2);
    FR(1:24,1) = 1/24;
    

    for t = 1:TMAX
        for L = 1:TMAX+24
            if t-L > 0
                FR(t,2) = FR(t,2) + FR(t-L,1) * UHM(L);
            end
        end
    end
    
    for t = 1:UH_DAY
        UH_river(t) = sum(FR(24*t-23:24*t,2));
    end
        
end
%--------------------------------------------------------------------------


%--------------- HRU's UH represented by Gamma distribution ---------------
UH_HRU_direct = zeros(KE,1);
for i = 1:KE
    x = linspace(24*(i-1),24*i,1001);
    hruh_fun = 1/(1/K)/gamma(N)*(x/(1/K)).^(N-1).*exp(-x/(1/K));
    pinteg = sum(hruh_fun)*(x(2)-x(1));
    UH_HRU_direct(i) = pinteg;
end
UH_HRU_base    = zeros(KE,1);
UH_HRU_base(1) = 1;
%--------------------------------------------------------------------------


%-------- Combined UH for HRU's response at the watershed outlet ----------
UH_direct = zeros(1,KE+UH_DAY-1);
UH_base   = zeros(1,KE+UH_DAY-1);
for k = 1:KE
    for u = 1:UH_DAY
        UH_direct(k+u-1) = UH_direct(k+u-1) + UH_HRU_direct(k) * UH_river(u);
        UH_base(k+u-1)   = UH_base(k+u-1) + UH_HRU_base(k) * UH_river(u);
    end
end
UH_direct = UH_direct/sum(UH_direct);
UH_base   = UH_base/sum(UH_base);
%--------------------------------------------------------------------------


%------------ Make Convolution for watershed outlet total flow ------------
directflow = zeros(length(inflow_direct),1);
baseflow   = zeros(length(inflow_direct),1);
for i = 1:length(inflow_direct)
    for j = 1:KE+UH_DAY-1
        if i-j+1 >= 1
            directflow(i) = directflow(i) + UH_direct(j) * inflow_direct(i-j+1);
            baseflow(i)   = baseflow(i) + UH_base(j) * inflow_base(i-j+1);
        end
    end
end
runoff = directflow + baseflow;
%--------------------------------------------------------------------------





