function [syn] = boot(data,n,m)

% Bootstrap a dataset
%
% Inputs:
%   data: data set to be bootstrapped
%   n: number of bootstrap iterations
%   m: size of bootstrapped sample
% 
% Outputs:
%   syn: matrix of bootstrapped samples.

syn=zeros(m,n); % initialize output matrix

% Loop to bootstrap dataset n times
for i=1:n
    for j=1:m
        syn(j,i) = data(randi(length(data)));
    end
end

end

