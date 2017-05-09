function [Em,str] = modulusRenormPow2Fit(fitObj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Em0,str]=dataextract(); %data grabbed from Jason's excel file
a=fitObj.a;b=fitObj.b;c=fitObj.c;

Em=(Em0).*(a*(Em0).^(b)+c);

end

function [Em str]=dataextract()
x= [100	7.59E-03;
125.89	8.62E-03;
158.49	9.68E-03;
199.53	1.07E-02;
251.19	1.17E-02;
316.23	1.26E-02;
398.11	1.35E-02;
501.19	1.41E-02;
630.96	1.47E-02;
794.33	1.50E-02;
1000	1.51E-02;
1258.9	0.015;
1584.9	0.014689;
1995.3	0.014172;
2511.9	0.013465;
3162.3	0.012598;
3981.1	0.011604;
5011.9	0.010523;
6309.6	0.0094002;
7943.3	0.0082749;
10000	0.0071844];
Em=x(:,1)/1000;%convert to kPa
str=x(:,2);

end