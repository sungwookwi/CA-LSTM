function [snow_outflow, snow_swe, snow_melt] = snow_ddf1(Prcp, Tavg, SnowPar, IniState)
%#codegen
% Snow process model using a simple Degree Day Factor model
%
% Inputs
%   Prcp     : Time series of daily precipitation
%   Tavg     : Time series of daily average temperature in Celsius
%   SnowPar  : Snow model parameters
%   IniState : Initial condition of snow storage
% 
% Snow Model Parameters
%	par(1) :  Degree-Day Factor (DDF)
%	par(2) :  Temperature threshold for snowfall (Ts)
%	par(3) :  Temperature threshold for snow melt (Tm)
% 
% Outputs
%   outflow :  Total outflow from the snowpack
%	swe     :  State of snow storage tank (= state of SWE)
%	melt    :  The amount of snow melt by DDF snow routine
%  
% codegen codes
%	cfg = coder.config('mex');
%   cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%   codegen snow_ddf1 -args {coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),zeros(1,3),0} -config cfg -report
% 
% Release Notes
% 	Written by Sungwook Wi(sungwookwi@gmail.com) 
%	Last Updated 01/08/2014 
% 


%------------------------------ SNOW parameters ---------------------------
DDF    = SnowPar(1);  % Degree-Day Factor
Ts     = SnowPar(2);  % Temperature threshold for snow
Tm     = SnowPar(3);  % Temperature threshold for snow melt
%--------------------------------------------------------------------------


%---------------------- Initialize snow storage state ---------------------
inisto_snow   = IniState;    % Initial state of snow storage
%--------------------------------------------------------------------------


%--------------------------- Output Arrays --------------------------------
snow_outflow = zeros(length(Prcp),1);    % total outflow from snow storage
snow_melt    = zeros(length(Prcp),1);    % snow melt
snow_swe     = zeros(length(Prcp),1);    % state of snow storage (SWE)
%--------------------------------------------------------------------------


%--------------------Execute model for every time step---------------------
for t = 1:length(Prcp)
    
    if Tavg(t) > Ts  % No snow, if temp is lower than the threshold
        snow    = 0;
        rain = Prcp(t);
    else              % Otherwise, we have snow for the precip
        snow    = Prcp(t);
        rain = 0;
    end
    
    % Updating snow storage with snow
    swe = inisto_snow + snow;
    
    % Determining wheather snow melt occurs and amount of snow melt using DDF
    if Tavg(t) > Tm  % Snow melt occurs if temp is above the base temp(Tb)
        melt = min(DDF*(Tavg(t)-Tm),swe);
    else             % No snow melt, otherwise
        melt = 0;
    end
    
    % Second snow storage update with melt
    swe       = swe - melt;
    inisto_snow  = swe;
    
    snow_outflow(t) = melt + rain;
    snow_melt(t)    = melt;
    snow_swe(t)     = swe;
    
end
%--------------------------------------------------------------------------





















