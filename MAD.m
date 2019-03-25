function [goodData,badData]=MAD(data,n)

%Funtion to calculate the mean absolute deviation (MAD)

%Inputs: 
%   data: a vector of data
%   n: MAD threshold for being an outlier
%
%Outputs: Vectors of sorted good and bad data

goodData=[];
badData=[];
%Calculate the MADs
m=median(data);
d_data=abs(data-m);
d_m=median(d_data);
MAD=d_data/d_m;
%Sort the data
for i=1:length(data)
    if MAD(i)<=n
        goodData=[goodData data(i)];
    else
        badData=[badData data(i)];
    end
end
        