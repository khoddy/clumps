function [cleanData,badData] = peirce(data)

% Performs peirce outlier test on array of data.
% Based on Peirces eq: 
%
%                       |xi-xm|max
%                R =   -----------
%                          stdev
%
% Where R is taken from a table of values.
%
% Currently only handles sample sizes between 3 and 12.
%
% Input: an array of values to be tested
% Outputs: arrays of good and bad data

% Basic stats for calculation
x_m = mean(data);
n = length(data);
sig = std(data);
d_data = abs(data-x_m);
x=[];
y=[];

load peirceTable

if n < 3
    disp('Sample size too small')
elseif n > 12
    disp('Sample size too large. Extend table of Peirce"s R values.')
else
    R = P(n-2,:);
    for i=1:length(R)        
        if R(i)>0 && length(y)==(i-1)
            y=[];
            x=[];
            d_max=R(i)*sig;
            for j=1:length(d_data) % loop to sort good and bad data
                if d_data(j)>d_max
                    y = [y data(j)]; %Bad data
                else
                    x = [x data(j)]; %Good data
                end
            end                      
        else    
        end
    end
end    
    
badData = y;    
cleanData = x;    