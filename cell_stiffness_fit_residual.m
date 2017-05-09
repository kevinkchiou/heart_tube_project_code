function [result] = cell_stiffness_fit_residual(Ecs,strains,Esats)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Use Stephanie's values to fit
[ec ece]=expt_cellgel_vals(1);
Es=ec(1,:);V = ec(2,:); W = max(ece);% Stephanie scaled values and error weights

%Use Engler's fit values


R2=zeros(numel(Ecs),numel(strains),numel(Esats));
for i=1:numel(Ecs)
    Ec=Ecs(i);
    for j=1:numel(strains)
        eps=strains(j);
        for k=1:numel(Esats)
            Esat=Esats(k);
            val=3*eps*compute_eshelby_trace_out(Ec,0.4,0.4);
            fmod_linear=min([Es/Esat;ones(1,numel(Es))]);
            v=val.*fmod_linear;
            R2(i,j,k)=sum(((V-v).*W).^2);
        end
    end
end

[min3,idx3]=min(R2,[],3);
[min2,idx2]=min(min3,[],2);
[~,idx1]=min(min2);

imin=idx1; %first index
jmin=idx2(imin); %second index
kmin=idx3(imin,jmin); %third index


Ec_min=Ecs(imin);eps_min=strains(jmin);Esat_min=Esats(kmin);
result.cell_stiffness=Ec_min;
result.strain_factor=eps_min;
result.saturating_stiffness=Esat_min;

[ec ece]=expt_cellgel_vals(1);
Gother=linspace(0.1,45,200);
figure(200);
L1=compute_eshelby_trace_out(Ec_min,0.4,0.4,Gother);
fmod_linear_Gother=min([linspace(0,max(Gother),numel(Gother))/Esat_min;ones(1,numel(Gother))]);
set(plot(log10(Gother),3*eps_min*L1.*fmod_linear_Gother,'r-'),'LineWidth',2);hold on; %factor of three when trace is taken
errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
%title('Constant stress source, spherical Eshelby cell');
xlabel('log_{10} ( E / E^* )');
%legend('Fit strain out','Experiment','Location','NorthEast');
hold off;
end

function f1=compute_eshelby_trace_out(E,nu_c,nu_t,Es)
%could in principle compute non-shperical or whatever, but sticking with
%this for now
if nargin<4
    Es=[0.3 0.84 2.6 10 40];
end
f1=((1-nu_t)/(1+nu_t))./(1+2.*Es*(1-2*nu_c)./(E*(1+nu_t)));
end

%Stephanie's values
function [val err]=expt_cellgel_vals(bool,norm)
if nargin<2
    norm=1/0.05;
end
if bool==0
    val=0;err=0;
else
    val=[0.3 0.84 2.6 10 40;norm*[0.022342 .049 .03649 .016 .0077]];
    dev=norm*[.017 .01009 .009 .0078 .012];
    err=[dev;-dev];
    err(1,val(2,:)-dev<0)=val(2,val(2,:)-dev<0); %prevent error bars < 0
end
end
