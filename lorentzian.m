function out = lorentzian(p,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if length(p)==3
    x0=p(1);g=p(2);scale=p(3);
elseif length(p)==2
    x0=p(1);g=p(2);
else
    error('Error: parameters are incorrect!\n');
end
out=scale/pi*(g./((x-x0).^2+g^2));
end