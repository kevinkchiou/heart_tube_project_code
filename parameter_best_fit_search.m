function [pars_out_3d pars_out_1d data3d data1d] = parameter_best_fit_search(figs,modeltype,alpha_range,gamma_range)
%PARAMETER_BEST_FIT_SEARCH() takes results from heart tube solutions and
%finds best fit results to experimental data. Various measures are used
%   Detailed fuck me goes here

E=[0.2 0.5 0.7 1 1.5 2.5]; %the values from experiment
if nargin<1
    figs=[1 2];
elseif length(figs)~=2
    figs(2)=figs(1)+1;
end
if nargin<2
    modeltype=0;
end
if modeltype==1
    E=E(1:end-1); %cut out high-stiffness cutoff
elseif modeltype==2
    E=E(2:end-1);
end
if nargin<3
    alpha = logspace(-1,1);
elseif length(alpha_range)~=2
    warning('Alpha:IncorrectSize','Wrong format for alpha_range argument!\n');
    warning('Alpha:IncorrectSize','Reverting to default values\n');
    alpha = logspace(-3,1);
else
    alpha = logspace(log10(alpha_range(1)),log10(alpha_range(2)),100);
end
if nargin<4
    gamma = linspace(0.001,0.10);
elseif length(gamma_range)~=2
    warning('Gamma:IncorrectSize','Wrong format for gamma_range argument!\n');
    warning('Gamma:IncorrectSize','Reverting to default values\n');
    gamma = linspace(0.001,0.07);
else
    gamma = linspace(gamma_range(1),gamma_range(2),100);
end

%heat map dataset for fit metric
data_hm_3d=zeros(length(alpha),length(gamma));
data_hm_1d=zeros(length(alpha),length(gamma));
for i=1:length(alpha);
    for j=1:length(gamma);
        %experimental data option set to zero since we are evaluating at
        %those points exclusively
        [v v1d]=threshold_wave_3d_solution(0,0,E,alpha(i),gamma(j));
        %now compute fit measure
        data_hm_3d(i,j)=compute_vel_fit_measure([E',v],modeltype);
        data_hm_1d(i,j)=compute_vel_fit_measure([E',v1d],modeltype);
    end
end
figure(figs(1));surf(gamma,log10(alpha),data_hm_3d);
figure(figs(2));surf(gamma,log10(alpha),data_hm_1d);
data3d=data_hm_3d;data1d=data_hm_1d;
[~,gidx3d]=min(min(data3d));[~,gidx1d]=min(min(data1d));
[m3d,aidx3d]=min(data3d(:,gidx3d));[m1d,aidx1d]=min(data1d(:,gidx1d));
if(~isempty(data3d(data3d<m3d)))
    error('oops!');
end
if(~isempty(data1d(data1d<m1d)))
    error('oops!');
end
pars_out_3d=[alpha(aidx3d) gamma(gidx3d)];
pars_out_1d=[alpha(aidx1d) gamma(gidx1d)];
end
