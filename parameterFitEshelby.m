function [pars,vmodel,vexp,Rsq] = parameterFitEshelby(fig,modeltype,alpharange,gammarange)
%PARAMETERFITESHELBY() takes results from heart tube solutions and
%finds best fit results to experimental data. Newer version that uses a
%nonlinear least squares approach instead of direct search

%single cell fit global variables
global Q;
global Es;
global Ec;
global eps_trans;

%single cell fit values - comsole (Jason) or fullspace (Kevin)
%Q=1;Es=0.9358;Ec=5.691;eps_trans=0.1575;%Jason's values
Q=0.1;Es=0.7323;Ec=0.7;eps_trans=1.05;%Kevin's values


E=sort([linspace(0.1,3,100),0.2,0.5,0.7,1,1.5,2.5]); %add values from experiment
Etemp=[0.2,0.5,0.7,1,1.5,2.5]; %trial E. only experimental values
if nargin<1
    fig=50;
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
    alpharange = 10.^[-1 1];
elseif length(alpharange)~=2
    warning('Alpha:IncorrectSize','Wrong format for alpha_range argument!\n');
    warning('Alpha:IncorrectSize','Reverting to default values\n');
    alpharange = 10.^[-3,1];
end
if nargin<4
    gammarange = [0.001 0.5];
elseif length(gammarange)~=2
    warning('Gamma:IncorrectSize','Wrong format for gamma_range argument!\n');
    warning('Gamma:IncorrectSize','Reverting to default values\n');
    gammarange = [0.001,0.5];
end
a=alpharange;g=gammarange;

x0=[1.30,0.0017];options=optimset('TolFun',1e-8,'TolX',1e-8);
%x0(2) will be multiplied by 1e-3 in computation; this is so finite
%difference changes dx make sense as a vector.

%check for initial guess
threshold_wave_3d_eshelby(fig,Etemp,x0(1),x0(2),2,[Q,Es,Ec,eps_trans]);
%pars = lsqnonlin(@eshelbyFitFunc,x0,[a(1) g(1)],[a(2) g(2)]);
pars=lsqnonlin(@eshelbyFitFunc,x0,[a(1) g(1)],[a(2) g(2)],options); %try without bounds
%non-derivative options
%pars=patternsearch(@sumEshelbyFitFunc,x0);
%pars=fminsearch(@sumEshelbyFitFunc,x0);


afit=pars(1);gfit=pars(2);
%replace check plot with full plot with all values of E
modeltemp=threshold_wave_3d_eshelby(fig,E,afit,gfit,2,[Q,Es,Ec,eps_trans]);
expttemp=expt_vel_vals(2);idx=zeros(length(expttemp),1);
for i=1:length(expttemp)
    [~,idx(i)]=min(abs(E-expttemp(1,i)));
end
vmodel=modeltemp(idx);
vexp=expttemp(2,:);

Rsq=sum((vmodel-vexp').^2);
end

function lsqvec=eshelbyFitFunc(x)
%eshelby fit parameters
global Q;
global Es;
global Ec;
global eps_trans;


E=[0.2,0.5,0.7,1,1.5,2.5]; %values from experiment
alpha=x(1);gamma=x(2); %input parameters to fit
v=threshold_wave_3d_eshelby(0,E,alpha,gamma,2,[Q,Es,Ec,eps_trans]);

%now compute the least squares metric
[ve ee]=expt_vel_vals(2);
de=ee(2,:)+ee(1,:); %error bar size
didx=de==0;de(didx)=mean(de(~didx));%correct for weird error bar of last point
lsqvec=(v-ve(2,:)')./de'; %vector output for error bars
%lsqvec=v-ve(2,:)'; %vector output for no error bars



%following is in case we require a finer E resolution for fzero
%{
idx=zeros(length(ve),1);
for i=1:length(ve)
    [~,idx(i)]=min(abs(E-ve(1,i)));
end
%lsqvec=(v(idx)-ve(2,:)').*de'; %output a vector

%following does without the error bars
lsqvec=(v(idx)-ve(2,:)').^2; %output a vector
%}
end

function lsqsum=sumEshelbyFitFunc(x)
lsqsum=sum(eshelbyFitFunc(x));
end

function [val err]=expt_vel_vals(bool)
if bool==0
    val=0;err=0;
elseif bool==1
    val=[0.2 0.5 0.7 1 1.5 2.5;1.98 10.29 13.56 18.22 27.03 0];
    err=[1.98 2.57 0.5 0 5.6 0;2.77 2.57 0.5 0 5.6 0];
elseif bool==2 %not so bool anymore :(
    val=[0.2 0.5 0.7 1 1.5 2.5;1.98 10.29 13.56 18.22 27.03 0];
    err=[1.98 2.57 0.5 0 5.6 0;2.77 2.57 0.5 0 5.6 0];
    %val=[0.2 0.5 0.7 1 1.5 2.5;1.98 0 13.56 18.22 27.03 0];
    %err=[1.98 0 0.5 0 5.6 0;2.77 0 0.5 0 5.6 0];
end
end
