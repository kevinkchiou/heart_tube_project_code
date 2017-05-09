function [I_1 I_2 I_3] = Tr_s2_factors(dr,dt,E,nu,Gamma,eta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Dpar=E*(1-nu) / (Gamma*(1+nu)*(1-2*nu));
Dperp=E/(2*Gamma*(1+nu));tau = E/(2*eta*(1+nu));

%I_i = the prefactor to v_i
%I_ji = the 1/q^j integral prefactor to v_i
result = 2*E^2/(3*((1+nu)*Gamma)^2)*I_factors(dr,dt,Dpar,0,0)+...
    E^2/(2*((1+nu)*Gamma)^2)*I_factors(dr,dt,Dperp,tau,0);
I_01 = result(1);I_02 = result(2);I_03 = result(3);
result = E^2/((1+nu)^2*Gamma*eta)*I_factors(dr,dt,Dperp,tau,2);
I_21 = result(1);I_22 = result(2);I_23 = result(3);
result=1/2*(E/((1+nu)*eta))^2 * I_factors(dr,dt,Dperp,tau,4);
I_41 = result(1);I_42 = result(2);I_43 = result(3);

%can't do the following 3 lines straight up - they have different
%contractions with Q. Need to separate and combine at higher up.
I_1 = I_01+I_21+I_41;
I_2 = I_02+I_22+I_42;
I_3 = I_03+I_23+I_43;

end

function [I1,I2,I3] = I_factors(r,t,D,tau,type)
%full computations for analytical form done via mathematica notebook
D=2*D; tau = tau/2; %since it is s^2 get extra factors of 2
a = D.*t./r.^2; sa = sqrt(a);ai = r.^2 / (4*D.*t); sai = sqrt(ai);
expa=exp(-ai/4);erfa=erf(sai/2);
if type==4
    I1=(expa.*(sai.^3+20*sai + 210*sa) - ...
        15*sqrt(pi).*(14*a-1).*erfa)./(8*pi^(3/2).*r.^3);
    I2=(-expa.*(2.*sai + 30*sa) + ...
        3*sqrt(pi)*(10*a - 1).*erfa)./(8*pi^(3/2).*r.^3);
    I3=(6*expa.*sa - ...
        (6*a-1)*sqrt(pi).*erfa)./(8*pi^(3/2).*r.^3);
elseif type==2
    I1=(-expa.*(sai.^7+14*sai.^5 + 140*sai.^3 + 840*sai) +... 
        840*sqrt(pi).*erfa)./(32*pi^(3/2).*r.^5);
    I2=(expa.*(sai.^5+10*sai.^3+60*sai) - ...
        60*sqrt(pi)*erfa)./(16*pi^(3/2).*r.^5);
    I3=(-expa.*(sai.^3+6*sai ) + ...
        6*sqrt(pi).*erfa)./(8*pi^(3/2).*r.^3);
elseif type==0
    I1=(expa).*sai.^11./(128*pi^(3/2)*r.^7);
    I2=(-expa).*sai.^9./(64*pi^(3/2)*r.^7);
    I3=(expa).*sai.^7./(32*pi^(3/2)*r.^7);
else
    fprintf('Selected type = %f!\n',type);
    fprintf('Type must = {0,2,4}, depending on q-integral!\n');
    error('Incorrect type!\n');
    %the type is the power of q in the denominator in the original 3d
    %Fourier transform ie in Int[d^3q / (2pi)^2 * 1/q^n exp(-iqx - Dq^2t)]
    %it should correspond to type = n
end
if tau~=0
    I1=I1*exp(-t/tau);
    I2=I2*exp(-t/tau);
    I3=I3*exp(-t/tau);
end
end