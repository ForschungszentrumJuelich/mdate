function [data] = getFOVdata(XC, TY, Ysum, X, T, amplitude)

%% if we choose nearest elements of X and T in XC and TY
[isInXC, locXC] = ismembertol(X, XC);
data = Ysum(locXC,:);
    
if exist('amplitude', 'var')
    data = data * amplitude / gmax(data);
end

[isInT, locTY] = ismembertol(T, TY');
data = data(:,locTY);


    
%% if we interpolate:

% interpolationMethod = 'linear';
% % in x
% data = interp1(XC, Ysum, X, interpolationMethod);
% 
% 
% % in t
% data = interp1(TY, data', T, interpolationMethod)';
% 
% if find(isnan(data))
%     a=1;
% end

% 
% normFactor =  amplitude / gmax(data);
% data = data * normFactor;

end