function [pars dat] = parameter_search_eshelby(figs,modeltype,alpharange,gammarange)
%PARAMETER_BEST_FIT_SEARCH() takes results from heart tube solutions and
%finds best fit results to experimental data. Various measures are used
%   Detailed fuck me goes here

E=[0.2 0.5 0.7 1 1.5 2.5]; %the values from experiment
if nargin<1
    figs=50;
end
if nargin<2
    modeltype=0;
end
if modeltype==1
    E=E(1:end-1); %cut out high-stiffness cutoff
elseif modeltype==2
    E=E(2:end-1); %cut out high and low end velocity cutoffs
end
if nargin<3
    alpha = logspace(-1,1);
elseif length(alpharange)~=2
    warning('Alpha:IncorrectSize','Wrong format for alpha_range argument!\n');
    warning('Alpha:IncorrectSize','Reverting to default values\n');
    alpha = logspace(-3,1);
else
    alpha = logspace(log10(alpharange(1)),log10(alpharange(2)),100);
end
if nargin<4
    gamma = linspace(0.001,0.10);
elseif length(gammarange)~=2
    warning('Gamma:IncorrectSize','Wrong format for gamma_range argument!\n');
    warning('Gamma:IncorrectSize','Reverting to default values\n');
    gamma = linspace(0.001,0.07);
else
    gamma = linspace(gammarange(1),gammarange(2),100);
end

%heat map dataset for fit metric
dat=zeros(length(alpha),length(gamma));
for i=1:length(alpha);
    for j=1:length(gamma);
        %experimental data option set to zero since we are evaluating at
        %those points exclusively
        v=threshold_wave_3d_eshelby(0,E,alpha(i),gamma(j),2);
        %now compute fit measure
        dat(i,j)=compute_vel_fit_measure([E',v],modeltype);
    end
end
figure(figs(1));surf(gamma,log10(alpha),dat);
xlabel('\Gamma');ylabel('log_{10} \alpha');zlabel('Residual');
[~,gidx]=min(min(dat));
[m,aidx]=min(dat(:,gidx));
if(~isempty(dat(dat<m)))
    error('oops!');
end
pars=[alpha(aidx) gamma(gidx)];
end


