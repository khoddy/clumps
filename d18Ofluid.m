function [d18Of, d18OfERR] = d18Ofluid(d18O,T,Terr)

% Keith Hodson - 08232017
%This function takes the d18O (VPDB) and clumped isotope temperature and
%calculates the parent fluid d18O and error (VSMOW). Uses Kim and
%O'Neil calibration for temperature dependent fractionation.

%Currently propagates the 1SE error in temperature, would probably be
%better if it used the SE value for D47 instead...

%Convert carb d18O into VSMOW:

d18O_SMOW = (1.03091*d18O)+30.91;

%Get fractionation factor (alpha):
alpha = exp(((18.03*(1000/(T+273.15)))-32.42)/1000);

%Get fluid d18O:
d18Of = ((d18O_SMOW+1000)/alpha)-1000;

%Propogate error from Terr:
alpha_err = exp(((18.03*(1000/(T+Terr+273.15)))-32.42)/1000);
d18Of_err = ((d18O_SMOW+1000)/alpha_err)-1000;
d18OfERR = (d18Of - d18Of_err);