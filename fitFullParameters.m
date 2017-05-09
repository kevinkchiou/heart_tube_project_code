function [pars,vmodel,convRatio,vexp,Rsq] = fitFullParameters(fig,modeltype,alpharange,gammarange)
%PARAMETERFITESHELBY() takes results from heart tube solutions and
%finds best fit results to experimental data. Newer version that uses a
%nonlinear least squares approach instead of direct search

%single cell fit global variables
global Q;
global eps_trans;
global err_select;

%single cell fit values - comsole (Jason) or fullspace (Kevin)
%Q=1;Es=0.9358;Ec=5.691;eps_trans=0.1575;%Jason's values
Q=0.1;eps_trans=1.;%Es=1;Ec=1;%Kevin's values

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
    %alpharange = 10.^[-2 0];
    alpharange = 10.^[-1 1];
elseif length(alpharange)~=2
    warning('Alpha:IncorrectSize','Wrong format for alpha_range argument!\n');
    warning('Alpha:IncorrectSize','Reverting to default values\n');
    alpharange = 10.^[-3,1];
end
if nargin<4
    gammarange = [0.0005 0.05];
elseif length(gammarange)~=2
    warning('Gamma:IncorrectSize','Wrong format for gamma_range argument!\n');
    warning('Gamma:IncorrectSize','Reverting to default values\n');
    gammarange = [0.001,0.002];
end
a=alpharange;g=gammarange;s=[0.25 4.0];c=[0.25 4.0];

%fitting a number of quantities - x0 = [alpha,gamma,Ec,Es]. different
%starting points find different minima for lsqnonlin

x0=[0.6,0.0017,0.7,0.7];err_select=1; %try a different starting point

options=optimset('TolFun',1e-4,'TolX',1e-4);
%x0(2) will be multiplied by 1e-3 in computation; this is so finite
%difference changes dx make sense as a vector.

%check for initial guess
threshold_wave_3d_eshelby(fig,Etemp,x0(1),x0(2),2,[Q,x0(3),x0(4),eps_trans]);
%pars = lsqnonlin(@eshelbyFitFunc,x0,[a(1) g(1)],[a(2) g(2)]);
pars1=lsqnonlin(@eshelbyFitFunc1,x0,[a(1) g(1) s(1) c(1)],[a(2) g(2) s(2) c(2)],options); %try without bounds
pars2=lsqnonlin(@eshelbyFitFunc2,x0,[a(1) g(1) s(1) c(1)],[a(2) g(2) s(2) c(2)],options); %try without bounds
%pars3=lsqnonlin(@eshelbyFitFuncStrain1,x0,[a(1) g(1) s(1) c(1)],[a(2) g(2) s(2) c(2)],options);
%pars4=lsqnonlin(@eshelbyFitFuncStrain2,x0,[a(1) g(1) s(1) c(1)],[a(2) g(2) s(2) c(2)],options);

%non-derivative options
%pars=patternsearch(@sumEshelbyFitFunc,x0);
%pars=fminsearch(@sumEshelbyFitFunc,x0);


afit1=pars1(1);gfit1=pars1(2);Esfit1=pars1(3);Ecfit1=pars1(4);
afit2=pars2(1);gfit2=pars2(2);Esfit2=pars2(3);Ecfit2=pars2(4);
%afit3=pars3(1);gfit3=pars3(2);Esfit3=pars3(3);Ecfit3=pars3(4);
%afit4=pars4(1);gfit4=pars4(2);Esfit4=pars4(3);Ecfit4=pars4(4);
%replace check plot with full plot with all values of E
[modeltemp1 eps1 conv1]=threshold_wave_3d_eshelby(0,E,afit1,gfit1,1,[Q,Esfit1,Ecfit1,eps_trans]);
[modeltemp2 eps2 conv2]=threshold_wave_3d_eshelby(0,E,afit2,gfit2,2,[Q,Esfit2,Ecfit2,eps_trans]);
%[modeltemp3 eps3 conv3]=threshold_wave_3d_eshelby(0,E,afit3,gfit3,1,[Q,Esfit3,Ecfit3,eps_trans]);
%[modeltemp4 eps4 conv4]=threshold_wave_3d_eshelby(0,E,afit4,gfit4,2,[Q,Esfit4,Ecfit4,eps_trans]);
plotResults(fig,E,[modeltemp1,modeltemp2],[eps1,eps2]);
%plotResults(fig,E,[modeltemp1,modeltemp2,modeltemp3,modeltemp4],[eps1,eps2,eps3,eps4]);
expttemp=expt_vel_vals(2);idx=zeros(length(expttemp),1);
for i=1:length(expttemp)
    [~,idx(i)]=min(abs(E-expttemp(1,i)));
end
vmodel1=modeltemp1(idx);
vmodel2=modeltemp2(idx);
%vmodel3=modeltemp3(idx);
%vmodel4=modeltemp4(idx);
vexp=expttemp(2,:);

%see if the two differ
[test11 test12]=threshold_wave_3d_eshelby(0,expttemp(1,:),afit1,gfit1,1,[Q,Esfit1,Ecfit1,eps_trans]);
[test21 test22]=threshold_wave_3d_eshelby(0,expttemp(1,:),afit2,gfit2,2,[Q,Esfit2,Ecfit2,eps_trans]);
plotResults(fig+1,expttemp(1,:),[test11,test21],[test12,test22]);

Rsq=zeros(1,2);
%Rsq=zeros(1,4);
Rsq(:,1)=sum((vmodel1-vexp').^2);
Rsq(:,2)=sum((vmodel2-vexp').^2);
%Rsq(:,3)=sum((vmodel3-vexp').^2);
%Rsq(:,4)=sum((vmodel4-vexp').^2);

pars=[pars1',pars2'];
vmodel=[vmodel1,vmodel2];
%pars=[pars1',pars2',pars3',pars4'];
%vmodel=[vmodel1,vmodel2,vmodel3,vmodel4];

%convergence output
convRatio=[conv1,conv2];
%convRatio=[conv1,conv2,conv3,conv4];
end

function out=eshelbyFitFunc1(x)
out=eshelbyFitFunc(x,1);
end

function out=eshelbyFitFunc2(x)
out=eshelbyFitFunc(x,2);
end

function lsqvec=eshelbyFitFunc(x,mod)
%eshelby fit parameters
global Q;
global eps_trans;
global err_select;

[ve ee]=expt_vel_vals(2); %ventricle data
%[ve ee]=expt_vel_vals(3); %additional upper points set to zero
%[ve ee]=expt_vel_vals(4); %additional upper and lowest points set to zero

%E=ve(1,:); %values from experiment
E=sort([linspace(0.1,3,20),ve(1,:)]); %slightly more values than from just experiment
alpha=x(1);gamma=x(2);Es=x(3);Ec=x(4); %input parameters to fit
vtemp=threshold_wave_3d_eshelby(0,E,alpha,gamma,mod,[Q,Es,Ec,eps_trans]);

Eidx=zeros(size(ve(1,:)));
for i=1:numel(Eidx)
    [~,Eidx(i)]=min(abs(E(:)-ve(1,i)));
end
v=vtemp(Eidx(:)); %pick out relevant velocity values

%now compute the least squares metric
if err_select==1
    de=ee(2,:)+ee(1,:); %error bar size
    didx=de==0;de(didx)=mean(de(~didx));%correct for error bars with small variance
    de(end)=0.0001*mean(de(~didx)); %make last point more important;
    lsqvec=(v-ve(2,:)')./de'; %vector output for error bars
elseif err_select==0
    lsqvec=v-ve(2,:)'; %vector output for no error bars
else
    error('fuck me');
end


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

function lsq=eshelbyFitFuncStrain1(x)
lsq=eshelbyFitFuncStrain(x,1);%second argument is saturating stress or not
end

function lsq=eshelbyFitFuncStrain2(x)
lsq=eshelbyFitFuncStrain(x,2);%second argument is saturating stress or not
end

function lsqvec=eshelbyFitFuncStrain(x,mod)
global Q;
global eps_trans;
global err_select;

E=[0.2,0.7,1,1.5,2.5]; %values from experiment
alpha=x(1);gamma=x(2);Es=x(3);Ec=x(4); %input parameters to fit
%use strain fit instead of standard fit
[~,v]=threshold_wave_3d_eshelby(0,E,alpha,gamma,mod,[Q,Es,Ec,eps_trans]);

%now compute the least squares metric
[ve ee]=expt_strain_vals(3); %eliminate second term for ventricle data
if err_select==1
    de=ee(2,:)+ee(1,:); %error bar size
    didx=de==0;de(didx)=mean(de(~didx));%correct for weird error bar of last point
    de(end)=0.1*mean(de(~didx)); %make last point more important;
    lsqvec=(v-ve(2,:)')./de'; %vector output for error bars
else
    lsqvec=v-ve(2,:)'; %vector output for no error bars
end
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
elseif bool==3 %not so bool anymore :(
    val=[0.2 0.5 0.7 1 1.5 2.5 3.0 3.5 4.0 4.5 5.0 5.5;1.98 10.29 13.56 18.22 27.03 0 0 0 0 0 0 0];
    err=[1.98 2.57 0.5 0 5.6 0 1 1 1 1 1 1;2.77 2.57 0.5 0 5.6 0 1 1 1 1 1 1];
elseif bool==4 %not so bool anymore :(
    val=[0.2 0.5 0.7 1 1.5 2.5 3.0 3.5 4.0 4.5 5.0 5.5;0 0 13.56 18.22 27.03 0 0 0 0 0 0 0];
    err=[1.0 2.57 0.5 0 5.6 0 1 1 1 1 1 1;1.0 2.57 0.5 0 5.6 0 1 1 1 1 1 1];
end
end

function [val err]=expt_strain_vals(bool)
if bool==0
    val=0;err=0;
elseif bool==1
    %val=[0.2 0.5 0.7 1 1.5 2.5;0.125*[0.27 0.2 0.72 1.0 0.8 0.5]];
    %err=.125*[.15 .12 .22 .03 .04 .04;.14 .13 .25 .06 .05 .05];
    val=[0.2 0.5 0.7 1 1.5 2.5;0.125*[0.27 0.2 0.72 1.0 0.8 0]];
    err=.125*[.15 .12 .22 .03 .04 0;.14 .13 .25 .06 .05 0];
elseif bool==2
    val=[0.2 0.5 0.7 1 1.5 2.5;0.0287 0 .1134 0.1104 0.0614 0];
    err=[0.0287 0 0.0758 0.0203 0.0039 0;0.0356 0 0.0758 0.0203 0.0039 0];
    %val=[0.2 0.5 0.7 1 1.5 2.5;0 0 .1134 0.1104 0.0614 0];
    %err=[0 0 0.0758 0.0203 0.0039 0;0 0 0.0758 0.0203 0.0039 0];
elseif bool==3
    val=[0.2 0.7 1 1.5 2.5;0.0287 .1134 0.1104 0.0614 0];
    err=[0.0287 0.0758 0.0203 0.0039 0;0.0356 0.0758 0.0203 0.0039 0];
    %val=[0.2 0.5 0.7 1 1.5 2.5;0 0 .1134 0.1104 0.0614 0];
    %err=[0 0 0.0758 0.0203 0.0039 0;0 0 0.0758 0.0203 0.0039 0];
end
end

function [val err]=expt_cellgel_vals(bool,norm)
if nargin<2
    norm=1/0.05;
end
if bool==0
    val=0;err=0;
elseif bool==1
    val=[0.3 0.84 2.6 10 40;norm*[0.022342 .049 .03649 .016 .0077]];
    dev=norm*[.017 .01009 .009 .0078 .012];
    err=[dev;-dev];
    err(1,end)=-val(2,end); %adjusts lower bar of last entry
end
end

function out=plotResults(fig,E,v,eps)
%[ec ece]=expt_cellgel_vals(1);
[es ese]=expt_strain_vals(2);
[ev eve]=expt_vel_vals(2);

szv=size(v);sze=size(eps);
if szv(2)>szv(1)
    v=v';szv=size(v);
end
if sze(2)>szv(1)
    eps=eps';sze=size(eps);
end
%error check
if sze(1)~=length(E) || szv(1)~=length(E)
   error('Error: incorrect data sizes!'); 
end
if szv(2)==2
    cmp(1,:)=[1 0 0];cmp(2,:)=[0 0 1];
end
if szv(2)==4
    cmp(1,:)=[1 0 0];cmp(2,:)=[0 0 1];cmp(3,:)=[0 1 0];cmp(4,:)=[0 0 0];
end
figure(fig);close(fig);set(figure(fig),'Position',[0 0 600 1000]);
for i=1:szv(2)
    vel=v(:,i);strain=eps(:,i);subplot(2,1,1);%plot(G,vel,'r.');
    hold on;plot(E,vel,'Color',cmp(i,:));
    errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'bs');
    axis([0 3 0 40]);
    xlabel('E / E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
    subplot(2,1,2);hold on;plot(E,strain,'Color',cmp(i,:));
    errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'bs');
    xlabel('E / E^*');ylabel('\epsilon (% strain)');
    subplot(2,1,1);
end
out=fig;
end