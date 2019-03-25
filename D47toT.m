function [T, Terr] = D47toT(D47,Cal,err)
% Turns a D47 measurement into temperature in C using the calibration 
% of your choice. Also calculates a temperature error from D47 errors.
% If error not known, input 0 for err.

% NOTE: Be sure your acid fractionation offsets are appropriate for your 
% chosen calibration!

% 'K' = Kelson et al. (2017) 90C *** Probably the one you want.
% 'k90' = Kluge (2015) 90C
% 'k25' = Kluge (2015) 25C
% 'z' = Zaarur et al. (2013)
% 'ph' = Passey and Henkes (2012) - disabled...
% 'g' = Ghosh et al. (2006)
% 'ds' = Dennis and Schrag (2010)

%=========================================================================

% Kluge et al 2015 (90C)
if strcmp(Cal,'k90') == 1
    T = sqrt((0.405*(10^5))/(D47-0.167));
    Terr = abs(sqrt((0.405*(10^5))/(D47+err-0.167))-T);
% Kluge et al 2015 (25C)
elseif strcmp(Cal,'k25') == 1    
    T = sqrt((0.397*(10^5))/(D47-0.167));
    Terr = abs(sqrt((0.397*(10^5))/(D47+err-0.248))-T);
% Zaarur et al. (2013)    
elseif strcmp(Cal,'z') == 1
    T = sqrt((0.555*(10^5))/(D47-0.078));
    Terr = abs(sqrt((0.555*(10^5))/(D47+err-0.078))-T);
% Passey and Henkes (2012) NEED TO ADD SOLVER FUNCTION!?!
%elseif strcmp(Cal,'ph') == 1
%D47 = -3.407*10^9/T^4 + 2.365*10^7/T^3 - 2.607*10^3/T^2 - 5.88/T + 0.28;
% Ghosh et al. (2006)
elseif strcmp(Cal,'g') == 1
    T = sqrt((0.636*(10^5))/(D47-0.005));
    Terr = abs(sqrt((0.636*(10^5))/(D47+err-0.005))-T);
% Dennis and Schrag (2010)
elseif strcmp(Cal, 'ds') == 1
    T = sqrt((0.018*(10^5))/(D47-0.292));
    Terr = abs(sqrt((0.018*(10^5))/(D47+err-0.292))-T);
% Kelson et al 2017 (90C)
elseif strcmp(Cal,'K') == 1
    T = sqrt((0.0417*(10^6))/(D47-0.139));
    Terr = abs(sqrt((0.0417*(10^6))/(D47+err-0.139))-T);
else
    sprintf('Calibration not recognized')
end
T = T-273.15;
end