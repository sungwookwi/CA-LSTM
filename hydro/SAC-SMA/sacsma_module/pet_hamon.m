function pet = pet_hamon(S_Date, E_Date, Tavg, Latitude, Coeff)
%#codegen
% Calculation of Potential Evapotranspiration based on Hamon method
%
% Inputs
%   S_Date    = Start Date of Simulation [year, month, day]
%   E_Date    = End Date of Simulation [year month day]
%   Tavg      = Time series of daily average temperature(C) for the specified simulation period
%   Latitude = Lough estimate of latitude of the area
%   Coeff     = Hamon PET adjustment factor
%
% codegen codes for generating MEX function
%   cfg = coder.config('mex');
%   cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%   codegen pet_hamon -args {zeros(1,3),zeros(1,3),coder.typeof(0,[200000 1],[1 0]),0,0} -config cfg -report
%
% written by Sungwook Wi
% sungwookwi@gmail.com
% 01/07/2014
%

% See References below for details
%   Lu Jianbiao, Ge Sun, Steven G. McNulty, Devendra M. Amatya (2005), A
%       comparison of six potential evaportranspiration methods for
%       regional use in the southeastern United States. Journal of the
%       American Water Resources Association, Vol. 41, No. 3., pp. 621-633,
%       doi:10.1111/j.1752-1688.2005.tb03759.x
%   Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahla,
%       Hsini Wua, Robert M. Schoolfieldb (1995), A model comparison for
%       daylength as a function of latitude and day of year Ecological
%       Modelling, Volume 80, Issue 1, Pages 87–95

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Daylight hour calculation based on CBM model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yrmat = (S_Date(1):E_Date(1))';
daylighthr_mat = zeros(length(Tavg),1);
s_ind = 1;
for iy = 1:length(yrmat)
    
    %------------------- Julian date generation ---------------------------
    if iy == 1
        if ( S_Date(2) <= 2 ) % January & February
            year  = S_Date(1) - 1.0;
            month = S_Date(2) + 12.0;
        else
            year = S_Date(1);
            month = S_Date(2);
        end
        day = S_Date(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year = S_Date(1) - 1.0;
        month = 1 + 12.0;
        day = 1;
        jd_2 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;  
        
        str_jdate = jd_1 - jd_2 + 1;
    else
        str_jdate = 1;
    end
    
    if iy == length(yrmat)
        if ( E_Date(2) <= 2 ) % January & February
            year  = E_Date(1) - 1.0;
            month = E_Date(2) + 12.0;
        else
            year = E_Date(1);
            month = E_Date(2);
        end
        day = E_Date(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year = E_Date(1) - 1.0;
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
    %----------------------------------------------------------------------
    
    var_theta   = 0.2163108 + 2*atan(0.9671396*tan(0.0086*(jdate-186)));
    var_pi      = asin(0.39795*cos(var_theta));
    daylighthr  = 24-24/pi*acos((sin(0.8333*pi/180)+sin(Latitude*pi/180)*sin(var_pi))./(cos(Latitude*pi/180)*cos(var_pi)));
    daylighthr_mat(s_ind:s_ind+length(jdate)-1) = daylighthr';
    
    s_ind = s_ind + length(jdate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Computing Potential Evapotranspiration based on Hamon euqation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
esat = 0.611*exp(17.27*Tavg./(237.3+Tavg));
pet  = Coeff*29.8*daylighthr_mat.*(esat./(Tavg+273.2));







