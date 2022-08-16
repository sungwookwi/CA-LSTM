function [snow_outflow, snow_melt, snow_swe, new_IniState, INtot] = snow_snow17_fix(S_date, E_date, Prcp, Tavg, Elev, Par, IniState)
%#codegen
% Snow process model based on SNOW17 (Anderson, 1976)
%
% Inputs
%   S_date   : Start date of simulation
%   E_date   : End date of simulation
%   Prcp     : Time series of daily precipitation
%   Tavg     : Time series of daily average temperature in Celsius
%   Elev     : Average elevation of the area
%   Par      : Snow model parameters
%   IniState : Initial condition of snow storage
%
% Snow Model Parameters
%	Par(1)   SCF    - Multiplying factor which adjusts Precipitation (accounts for gage snow catch deficiencies)
%	Par(2)   PXTEMP - Temperature that separates rain from snow, deg C
%	Par(3)   MFMAX  - Maximum melt factor during non-rain periods - assumed to occur on June 21
%	Par(4)   MFMIN  - Minimum melt factor during non-rain periods - assumed to occur on Dec 21
%	Par(5)   UADJ   - Average wind function during rain-on-snow periods
%	Par(6)   MBASE  - Base temperature for snowmelt computations during non-rain periods, deg C
%	Par(7)   TIPM   - Antecedent temperature index parameter (0.01 to 1.0)
%	Par(8)   PLWHC  - Percent liquid water holding capacity (maximum value allowed is 0.4)
%	Par(9)   NMF    - Maximum negative melt factor
%	Par(10)  DAYGM  - A constant daily rate of melt at the soil-snow interface
%	Par      TTI    - Temperature interval for mixture of snow and rain
%
% Outputs
%   outflow :  Total outflow from the snowpack
%	melt    :  The amount of snow melt by DDF snow routine
%	swe     :  State of snow storage tank (= state of SWE)
%
% Initial states
%	IniState(1)  W_i     -  Accumulated water equivalent of the ice portion of the snow cover (mm)
%   IniState(2)  ATI     -  Antecedent Temperature Index, deg C
%   IniState(3)  W_q     -  Liquid water held by the snow (mm)
%   IniState(4)  Deficit -  Heat Deficit, also known as NEGHS, Negative Heat Storage
% 
% codegen codes
%	cfg = coder.config('mex');
%   cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%   codegen snow_snow17 -args {zeros(1,3),zeros(1,3),coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),0,zeros(1,10),zeros(1,4)} -config cfg -report
%
% Release Notes
%	Written by Sungwook Wi(sungwookwi@gmail.com)
%	Last Updated 01/08/2014
%

% References
%	Anderson, E. A. (1976), A point energy and mass balance model of a snow
%       cover, NOAA Tech. Rep. NWS 19, 150 pp., Natl. Oceanic and Atmos.
%	Admin., Silver Spring, Md. Anderson, E. A. (2006), Snow accumulation
%       and ablation model - SNOW17,
%       http://www.nws.noaa.gov/oh/hrl/nwsrfs/users_manual/part2/_pdf/22snow17.pdf
%   


%                           SNOW-17 PARAMETERS
SCF    = Par(1);  % Multiplying factor which adjusts Precipitation
PXTEMP = Par(2);  % Temperature that separates rain from snow
% TTI    = Par(3);  % Temperature interval for mixture of snow and rain
MFMAX  = Par(3);  % Maximum melt factor during non-rain periods
MFMIN  = Par(4);  % Minimum melt factor during non-rain periods
UADJ   = Par(5);  % Average wind function during rain-on-snow periods
MBASE  = Par(6);  % Base temperature for snowmelt computations during non-rain periods
TIPM   = Par(7);  % Antecedent temperature index parameter
PLWHC  = Par(8);  % Percent liquid water holding capacity
NMF    = Par(9);  % Maximum negative melt factor
DAYGM  = Par(10); % A constant daily rate of melt at the soil-snow interface


%                     INITIALIZE STATE OF SNOW STORAGE
W_i     = IniState(1);  % Accumulated water equivalent of the ice
ATI     = IniState(2);  % Antecedent Temperature Index
W_q     = IniState(3);  % Liquid water held by the snow
Deficit = IniState(4);  % Heat Deficit



%                             OUTPUT ARRAYS
snow_outflow =  zeros(length(Prcp),1);    % total outflow from snow storage
snow_melt    =  zeros(length(Prcp),1);    % snow melt
snow_swe     =  zeros(length(Prcp),1);    % state of snow storage (SWE)
INtot     =  zeros(length(Prcp),1);    


new_IniState = nan(4,1);

%                           JULIAN DATE ARRAY
yrmat     = (S_date(1):E_date(1));
s_ind = 1;
jdate_mat      = zeros(length(Prcp),1);
yrmat2         = zeros(length(Prcp),1);
for iy = 1:length(yrmat) 
    
    % Julian date generation
    if iy == 1 
 
        if ( S_date(2) <= 2 ) % January & February
            year  = S_date(1) - 1.0;
            month = S_date(2) + 12.0;
        else
            year = S_date(1);
            month = S_date(2);
        end
        day = S_date(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year = S_date(1) - 1.0;
        month = 1 + 12.0;
        day = 1;
        jd_2 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        str_jdate = jd_1 - jd_2 + 1;
    else
        str_jdate = 1;
    end
    
    if iy == length(yrmat)
        
        if ( E_date(2) <= 2 ) % January & February
            year  = E_date(1) - 1.0;
            month = E_date(2) + 12.0;
        else
            year = E_date(1);
            month = E_date(2);
        end
        day = E_date(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year = E_date(1) - 1.0;
        month = 1 + 12.0;
        day = 1;
        jd_2 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        end_jdate = jd_1 - jd_2 + 1;
    else
        if ((mod(yrmat(iy),4) == 0 && mod(yrmat(iy),100) ~= 0) || mod(yrmat(iy),400) == 0)
            end_jdate = 366;
        else
            end_jdate = 365;
        end
    end
    
    jdate       = (str_jdate:end_jdate);
    numday      = length(jdate);
    jdate_mat(s_ind:s_ind+numday-1) = jdate';
    yrmat2(s_ind:s_ind+numday-1) = repmat(yrmat(iy),end_jdate-str_jdate+1,1);
    s_ind = s_ind + numday;   
end
% Date matrix with Julian day
Juliandd  =  [yrmat2 jdate_mat];  % [yr juliandate]



%                   EXECUTE MODEL FOR EVERY TIME STEP
dtt = 24; % time interval of temperature data
dtp = 24; % time interval of prcipitation data
for t = 1:length(Prcp)
    
    Ta = Tavg(t);     % Air temperature at this time step (deg C)
    Pr = Prcp(t);     % Precipitation at this time step (mm)
    

    %-------------------- FORM OF PRECIPITATION ---------------------------
    
    if Ta <= PXTEMP   % Air temperature is cold enough for snow to occur
        SNOW = Pr; RAIN = 0;
    else              % Air temperature is warm enough for rain
        SNOW = 0;  RAIN = Pr;
    end
    
%     if Ta >= (PXTEMP+TTI)   % All rain, no snow
%         SNOW = 0; 
%         RAIN = Pr;
%     elseif Ta <= PXTEMP % All snow, no rain
%         SNOW = Pr;  
%         RAIN = 0;       
%     else  % Linear mixture of snow and rain in interval TTI
%         snowfrac = -1/TTI * (Ta-PXTEMP) + 1;
%         SNOW = Pr*snowfrac;  
%         RAIN = Pr*(1-snowfrac);
%     end
    %----------------------------------------------------------------------


    %------------------ ACCUMULATION OF THE SNOW COVER --------------------
                 
    Pn  = SNOW*SCF;   % Water equivalent of new snowfall (mm)
    W_i = W_i + Pn;   % Water equivalent of the ice portion of the snow cover (mm)
    %E   = 0;          % Excess liquid water in the snow cover
    %----------------------------------------------------------------------
    

    %----- ENERGY EXCHANGE AT SNOW/AIR SURFACE DURING NON-MELT PERIODS ----

    % Seasonal variation in the non-rain melt factor
    DAYN = Juliandd(t,2);      % Current julian date
    if ((mod(Juliandd(t,1),4) == 0 && mod(Juliandd(t,1),100) ~= 0) || mod(Juliandd(t,1),400) == 0) % Leap year
        days=366;
        N_Mar21=DAYN-81;   % Day of year since March 21 (leap)
    else   % Not a leap year
        days=365;
        N_Mar21=DAYN-80;
    end
    Sv = (0.5*sin((N_Mar21 * 2 * pi)/days)) + 0.5;        % Seasonal variation
    Av = 1.0;                                             % Seasonal variation adjustment, Av=1.0 when lat < 54N
    Mf = (dtt/6) * ((Sv * Av * (MFMAX - MFMIN)) + MFMIN); % Seasonally varying non-rain melt factor
    
    % New snow temperature and heat deficit from new snow
    if Ta < 0
        T_snow_new = Ta;
    else
        T_snow_new = 0;
    end
    
    % Change in the heat deficit due to new snowfall (mm), 80 cal/g: latent
    % heat of fusion, 0.5 cal/g/C: specific heat of ice
    delta_HD_snow = -(T_snow_new*Pn)/(80/0.5);   

    % Change in heat deficit due to a temperature gradient
    delta_HD_T = NMF * (dtp/6) * (Mf/MFMAX) * (ATI - T_snow_new);
    
    % Update ATI(Antecedent Temperature Index)
    if Pn > (1.5*dtp)
        ATI = T_snow_new;      %Antecedent temperature index
    else
        TIPM_dtt = 1.0 - ((1.0 - TIPM)^(dtt/6));
        ATI = ATI + TIPM_dtt * (Ta - ATI);
    end
    ATI = min(ATI,0);
    %----------------------------------------------------------------------
   
    
    %----------------------------- SNOW MELT ------------------------------
                             
    T_rain = max(Ta,0);   % Temperature of rain (deg C), Ta or 0C, whichever greater
    if RAIN > (0.25 * dtp)  % Rain-on-Snow Melt
        stefan = 6.12*(10^(-10));  % Stefan-Boltzman constant (mm/K/hr)
        e_sat  = 2.7489*(10^8)*exp((-4278.63/(Ta+242.792)));  % Saturated vapor pressure at Ta (mb)
        % P_atm: Atmospheric pressure (mb) where elevation is in HUNDREDS
        % of meters (this is incorrectly stated in the manual)
        P_atm  = 33.86*(29.9-(0.335*(Elev/100))+(0.00022*((Elev/100)^2.4)));   
        term1 = stefan * dtp * (((Ta+273)^4)-(273^4));
        term2 = 0.0125 * RAIN * T_rain;
        term3 = 8.5 * UADJ * (dtp/6) * ((0.9*e_sat - 6.11) + (0.00057*P_atm*Ta));
        Melt = term1 + term2 + term3;
        Melt = max(Melt,0);
    elseif RAIN <= (0.25 * dtp) && (Ta > MBASE) % Non-Rain Melt
        Melt = (Mf * (Ta - MBASE) * (dtp/dtt)) + (0.0125 * RAIN * T_rain);
        Melt = max(Melt,0);
    else
        Melt = 0;
    end
    %----------------------------------------------------------------------
   
    
    %------------------ RIPENNESS OF THE SNOW COVER -----------------------
                    
    % W_i : water equivalent of the ice portion of the snow cover
    % W_q : liquide water held by the snow
    % W_qx: liquid water storage capacity
    % Qw  : Amount of available water due to melt and rain
    
    Deficit = max(Deficit + delta_HD_snow + delta_HD_T, 0);   % Deficit = heat deficit (mm)

    
    if Melt < W_i
        
        W_i = W_i-Melt;
        Qw = Melt + RAIN;
        W_qx = PLWHC * W_i;
        
        if Deficit>(0.33*W_i) % limits of heat deficit
            Deficit=0.33*W_i;
        end
        
        if (Qw + W_q) > (Deficit + Deficit*PLWHC + W_qx)
            % THEN the snow is RIPE
            E = Qw + W_q - W_qx - Deficit-(Deficit*PLWHC);  % Excess liquid water (mm)
            W_i = W_i + Deficit;  % W_i increases because water refreezes as heat deficit is decreased
            W_q = W_qx;  % fills liquid water capacity
            Deficit = 0;
            
        elseif (Qw >= Deficit) && ((Qw + W_q) <= ((Deficit*(1+PLWHC)) + W_qx))
            % THEN the snow is NOT yet ripe, but ice is being melted
            E = 0;
            W_i = W_i + Deficit; % W_i increases because water refreezes as heat deficit is decreased 
            W_q = W_q + Qw - Deficit;
            Deficit = 0;
            
        else
            % THEN the snow is NOT yet ripe
            E = 0;
            W_i = W_i + Qw; % W_i increases because water refreezes as heat deficit is decreased
            Deficit = Deficit - Qw;
        end
        
    else % Melt >= W_i
        
        Melt=W_i+W_q;
        W_i=0;
        W_q=0;
        Qw = Melt + RAIN;
        E = Qw;
        
    end
    
    if Deficit == 0
        ATI = 0;
    end
    %----------------------------------------------------------------------
     
    
    %------------------------- CONSTANT RELEASE ---------------------------
                                  
    % Constant daily amount of melt which takes place at the snow-soil
    % interface whenever there is a snow cover
    if W_i > DAYGM
        gmwlos = (DAYGM/W_i)*W_q;
        gmslos = DAYGM;
        gmro = gmwlos + gmslos;
        W_i = W_i - gmslos;
        W_q = W_q - gmwlos;
        
        E = E + gmro;
        SWE = W_i + W_q;
        
    else
        gmro = W_i+W_q;
        W_i=0;
        W_q=0;
        
        E = E + gmro;
        SWE = 0;
    end
    %----------------------------------------------------------------------
       
    
    %-------------------------- OUTPUTS -----------------------------------
    snow_outflow(t) =  E;      
    snow_melt(t)    =  Melt; 
    snow_swe(t)     =  SWE; 
    INtot(t)        =  Pn+RAIN;

end
    
new_IniState(1) = W_i;
new_IniState(2) = ATI;
new_IniState(3) = W_q;
new_IniState(4) = Deficit;
    
    
    
    
    
    
    
    
    
    
    
    
    
    