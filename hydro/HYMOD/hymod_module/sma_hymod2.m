function [direct, base, actevap, state_soil] = sma_hymod2(Water, PET, SWE, HYMOD_Par, IniState)
%#codegen 
% Soil Moisture Accounting Module based on HYMOD
% Daily total runoff simulation performed by hydrologic model HYMOD
% 
% Example Syntax
%	[direct base] = HYMOD_UMASS(pr, tas, pet, swe, par,[0, 100])
%   
% Inputs
%	Water      =  Daily input of available water
%	PET        =  Daily input of Potential Evapotranspiration
%	SWE        =  Daily input of Snow Water Equivalent
%	HYMOD_Par  =  HYMOD model parameters 
%	IniState   =  Initial Conditions of soil & groundwater storages
%
% Parameters
%	HYMOD_Par(1)  =  Coefficient of Permanent Wilting Point (Kpwp)
%	HYMOD_Par(2)  =  Exponent applied to the ratio of SM to PWP (ETexp)
%	HYMOD_Par(3)  =  Maximum storage capacity of soil moisture accounting tank (cmax)
%	HYMOD_Par(4)  =  Shape parameter of distribution function of storage capacity (B -> bexp)
%	HYMOD_Par(5)  =  Split parameter for Direct and Base runoff (alpha)
%	HYMOD_Par(6)  =  Groundwater release coefficient (Ks)
%	HYMOD_Par(7)  =  Lower reservoir maximum capacity (Lmax)
%
% Initial States of Storages
%	IniState(1)    =  Initial condition of soil moisture storage
%	IniState(2)    =  Initial condition of groundwater storage
%
% Outputs
%	direct  =  Daily Direct Runoff
%	base    =  Daily Baseflow
% 
% Other Possible Outputs
%	state_soil =  State of soil moisture account tank
%	actevap    =  Actual evaporation
%	effrain    =  Effective rainfall amount
%
% codegen codes
%	cfg = coder.config('mex');
%   cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%   codegen sma_hymod2 -args {coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),zeros(1,7),zeros(1,2)} -config cfg -report
% 
% Release Notes
% 	Written by Sungwook Wi(sungwookwi@gmail.com)  
% 

% References for more details
%	Vrugt, J. A., C. J. F. ter Braak, H. V. Gupta, and B. A. Robinson
%       (2009), Equifinality of formal (DREAM) and informal (GLUE) Bayesian
%       approaches in hydrologic modeling, Stoch. Environ. Res. Risk.
%       Assess., 23, 1011–1026, doi: 10.1007/s00,477–008–0274–y
%	Kollat, J. B., P. M. Reed, and T. Wagener (2012), When are
%       multiobjective calibration trade-offs in hydrologic models
%       meaningful?, Water Resour. Res., 48, W03520,
%       doi:10.1029/2011WR011534.
% 

%------------------------------ HYMOD parameters --------------------------
Kpwp   = HYMOD_Par(1);
ETexp  = HYMOD_Par(2);
cmax   = HYMOD_Par(3);
if HYMOD_Par(4) == 2  % Convert from scaled B (0-2) to unscaled b (0 - Inf)
    bexp = 10^6;
else
    bexp = log(1-(HYMOD_Par(4)/2))/log(0.5);
end
alpha  = HYMOD_Par(5);
Ks     = HYMOD_Par(6);
Lmax   = HYMOD_Par(7);
smax   = cmax/(1+bexp); % Maximum storage capacity
PWP    = Kpwp*smax; % Permanent Wilting Point
%--------------------------------------------------------------------------


%------------------------- Initial storage states -------------------------
soilwater    = IniState(1);    % Initial state of soil moisture storage
groundwater  = IniState(2);    % Initial state of groundwater storage

if soilwater > smax
    soilwater = smax;
end

if groundwater > Lmax
    groundwater = Lmax;
end
%--------------------------------------------------------------------------


%----------------------------- Output Arrays ------------------------------
simper     = length(Water);       % simulation period in days
state_soil = zeros(simper,1);    % state of soil moisture accounting tank
actevap    = zeros(simper,1);    % simulated actual evaporation
effrain    = zeros(simper,1);    % effective rainfall
%--------------------------------------------------------------------------


%--------------------Execute model for every time step---------------------
direct = nan(simper,1);  % simulated direct runoff
base   = nan(simper,1);  % simulated baseflow
for t = 1:simper  
    
    % Total amount of water available at a time step
    Wval        = Water(t);
    % Value of potential evaporation
    if SWE(t) > 0
        PETval      = 0.0; 
    else
        PETval      = PET(t);        
    end
    % State of soil moisture storage at the beginning of the time step 
    state_beg   = soilwater;            
 
    % Compute excess precipitation 
    cs_beg     = cmax*(1-(1-state_beg/smax)^(1/(bexp+1)));        % Critical storage capacity based on the beginning state
    OV1        = max(0,Wval+cs_beg-cmax);                         % Surface runoff
    Wval_inf   = Wval - OV1;                                      % Infiltrated water
    cs_end     = min(cs_beg+Wval_inf,cmax);                       % Updated critical storage capacity for current time step
    state_end  = smax*(1-(1-cs_end/cmax)^(1+bexp));               % State of soil moisture account by the end of the time step for precipitation event
    OV2        = max(Wval_inf-(state_end-state_beg),0);           % Surface, Subsurface and groundwater runoff based on the probability-distributed principle
    
    % Compute evaporation
    if state_end >= PWP
        Kpet = 1;
    else
        Kpet = (state_end/PWP)^ETexp;
    end
    et        = Kpet*PETval;   
    aet       = min(state_end,et);                           % Compute evaporation
    soilwater = state_end - aet;                             % Update state of soil moisture storage
    
    % Partition OV1 and OV2 into quick and slow flow component
    directflow = OV1 + alpha * OV2;  % Effective rainfall which comprises direct runoff
    effloss    = (1-alpha) * OV2;    % Effective abstraction (loss) contributing to groundwater storage
    
    % Change because of new parameter Lmax; when groundwater storage full, no groundwater recharge
    if (effloss+groundwater >= Lmax)
        directflow = directflow+(effloss+groundwater-Lmax);
        effloss = Lmax-groundwater;
    end
    
    % Groundwater storage update & groundwater release
    update_gwater  = (1-Ks)*groundwater + (1-Ks)*effloss;
    baseflow       = Ks/(1-Ks)*update_gwater;
    groundwater    = update_gwater; 
    
	% Intermediate states 
    state_soil(t) = soilwater+groundwater;         % State of soil moisture account at the beginning of time step
    actevap(t)    = aet;              % Computed actual evaporation
    effrain(t)    = OV1 + alpha*OV2;   % Effective rainfall derived by probability-distributed principle
    
    % Total outflow for the time step
    direct(t) = directflow;
    base(t)   = baseflow;
        
end


