function [ output_args ] = create_manuscript_plots(fig,fileoption)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<1 || numel(fig)==0
    fig=[1,2];
elseif numel(fig)==1
    fig=[fig,fig+1];
elseif numel(fig)==2
    fig=[fig(1),fig(2),fig(2)+1];
elseif numel(fig)>4
    fig=[fig(1),fig(2),fig(3),fig(4)];
end

E0=1.3;nu=0.4;Q=[1 0 0;0 -nu 0;0 0 -nu];
TrQ=trace(Q);
%TrQxy=Q(1,1)+Q(2,2);

r1=0.5318;Gamma1=0.0031;Q1=0.10;Esat1=0.700;Ec1=0.58;eps_trans1=1;%constant stress source (eshelby)
r2=0.5377;Gamma2=0.0029;Q2=0.10;Esat2=0.6823;Ec2=0.5947;eps_trans2=1;%constant stress source (eshelby)

if nargin<2;
    G=linspace(0.01,3,200);
    [v1 eps1]=threshold_wave_3d_eshelby(0,G,r1,Gamma1,1,[Q1,Esat1,Ec1,eps_trans1]);
    [v2 eps2]=threshold_wave_3d_eshelby(0,G,r2,Gamma2,2,[Q2,Esat2,Ec2,eps_trans2]);
    save('manuscript_plot_variables.mat','G','v1','eps1','v2','eps2');
elseif fileoption==1
    load('manuscript_plot_variables.mat');
elseif ischar(fileoption)
    load(fileoption);
end

%functions returning experimental values for errorbar() go here.
[es ese]=expt_strain_vals(3);
[ev eve]=expt_vel_vals(2);
[ec ece]=expt_cellgel_vals(1);
[ep epe]=expt_propagation_vals(1);

%plot tissue velocites and strains for each model
h=figure(fig(1));close(h);set(figure(fig(1)),'Position',[0 0 600 1000]);
subplot(2,1,1);plot(G,v1,'b-');hold on;plot(G,v2,'r-');
errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'ko');
xlabel('E / E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
subplot(2,1,2);plot(G,eps1,'b-');hold on;plot(G,eps2,'r-');
errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'ko');
xlabel('E / E^*');ylabel('Tissue strain ( Tr \epsilon )');hold off;
legend('Constant stress CM','Saturating stress CM','Experiment','Location','NorthEast');
legend boxoff;
%legend('Model','Experiment','Location','NorthEast');

%plot Eshelby out for each model (constant / saturating stress)
h=figure(fig(2));close(h);set(figure(fig(2)),'Position',[0 0 600 1000]);

subplot(2,1,1);
plot(G,2*Q1*ones(size(G)),'b:');hold on;
plot(G,2*Q1*min(1,G/Esat1),'r:');
legend('CE model','SE model','Location','NorthEast');
legend boxoff;
xlabel('E / E^*');ylabel('CM eigenstrain');
axis([0 3 0 0.3]);hold off;

Gother=logspace(-1,2,100);
%for plotting constant stress eshelby out
subplot(2,1,2);

%data=csvread('comsole_2d_3d_cellgel_jason.csv');
data=csvread('comsole_E_t_vs_s_sat_const.csv',1);
%semilogx(Gother,epsTrace2.*fmod_linear_Gother,'r-'); %factor of three when trace is taken
semilogx(data(:,1),data(:,3),'b-',data(:,1),data(:,2),'r-');hold on;
errorbar(ec(1,:),0.05*ec(2,:),0.05*ece(1,:),0.05*ece(2,:),'ko');
xlabel('E_{gel} ( kPa )');ylabel('Gel strain ( Tr \epsilon)');
axis([0.1,100,0,0.12]);hold off;
legend('FE simulation CE','FE simulation SE','Cell-on-gel experiment','Location','NorthEast');
legend boxoff;
%semilogx(Gother,2*Q1*ones(1,numel(Gother)),'b-',Gother,2*Q2*min([ones(1,numel(Gother));Gother/Esat2/E0]),'r-');
%legend('Constant stress','Saturating stress','Location','NorthEast');
%legend boxoff;
%axis([0.1 100 0 max(2*Q1,2*Q2)*1.5]);
%xlabel('E_{gel} ( kPa ) ');ylabel('CM eigenstrain');



h=figure(fig(3));close(h);set(figure(fig(3)),'Position',[0 0 600 1000]);
subplot(2,1,1);

[f1 f2]=compute_eshelby_factors(Gother/Ec1/E0,nu,nu);
epsTrace1=traceSelect(f1,f2,2*Q1*eps_trans1*Q,'partial');%epsTrace=traceSelect(f1,f2,2*Q1*eps_trans1,'full'); %functions for taking trace
%semilogx(Gother,epsTrace1,'b-');hold on; %factor of three when trace is taken
eps3dTrace1in=eshelby_mismatch_in(G/Ec1/E0,nu,nu);
[f1 ~]=compute_eshelby_factors(G/Ec1/E0,nu,nu);
eps3dTrace1out=3*f1;
%for plotting saturating stress eshelby out
[f1 f2]=compute_eshelby_factors(Gother/Ec2/E0,nu,nu);
epsTrace2=traceSelect(f1,f2,2*Q2*eps_trans2*Q,'partial');%epsTrace=traceSelect(f1,f2,2*Q2*eps_trans2,'full'); %functions for taking trace
[f1 ~]=compute_eshelby_factors(G/Ec2/E0,nu,nu);
eps3dTrace2out=3*f1;
eps3dTrace2in=eshelby_mismatch_in(G/Ec2/E0,nu,nu);
%following line is linear dependence of cell stress activation on E
fmod_linear_Gother=min([Gother/Esat2/E0;ones(1,numel(Gother))]);hold on;
%plotyy([Gother',Gother'],[3*TrQ*L1',(3*TrQ*L2.*fmod_linear_Gother)'],...
    %[Gother',Gother'],[2*Q1*ones(numel(Gother),1),2*Q2*min([ones(1,numel(Gother));Gother/Esat2])'],'semilogx','semilogx','semilogx','semilogx');


%plot(G,eps3dTrace1in,'r.',G,eps3dTrace1out,'b.',G,eps3dTrace1in.*eps3dTrace1out,'g.');
plot(G,eps3dTrace1in,'r-',G,eps3dTrace1out,'b-');
%semilogx(Gother,eps3dTrace1in,'r.',Gother,eps3dTrace1out,'b.');
legend('T^{in}','T^{out}','T^{in}*T^{out}','Location','NorthWest');
legend boxoff;
xlabel('E / E^*');
subplot(2,1,2);
set(errorbar(ep(1,:),ep(2,:),epe(1,:),epe(2,:),'^'),'Color',[40 150 100]/255);axis([0 3 0 1.1]);
legend('Propagation Probability','Location','NorthEast');
legend boxoff;


figure(fig(4));close(figure(fig(4)));set(figure(fig(4)),'Position',[0 0 400 400]);
Gshort=G(G<=1.5);
plot(Gshort,2*Q1*ones(size(Gshort)),'b:');hold on;
plot(Gshort,2*Q1*min(1,Gshort/Esat1),'r:');
legend('CE model','SE model','Location','NorthEast');
legend boxoff;
xlabel('E / E^*');ylabel('CM eigenstrain');
axis([0 1.5 0 0.25]);hold off;

plot2svg('./svgs/vel_strain.svg',figure(fig(1)));
plot2svg('./svgs/stress_strain.svg',figure(fig(2)));
plot2svg('./svgs/Tin_Tout_Prob.svg',figure(fig(3)));
plot2svg('./svgs/CE_SE_models.svg',figure(fig(4)));
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

function [val err]=expt_propagation_vals(bool)
if bool==0
    val=0;err=0;
elseif bool==1
    val=[0.2 0.5 0.7 1 1.5 2.5;0.2011 0.5 0.76 0.80 0.82 0];
    err=[0.16 0.1 0.3 0.26 0.05 0;0.16 0.1 0.3 0.26 0.05 0];
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

function out=traceSelect(f1,f2,Q,sel)
if nargin<4
    error('fuck me!');
end
trQ=trace(Q);
if strcmp(sel,'partial') || strcmp(sel,'Partial') || strcmp(sel,'PARTIAL')
    out=2*(f1*trQ + f2*(Q(1,1)+Q(2,2)-2/3*trQ));
elseif strcmp(sel,'full') || strcmp(sel,'Full') || strcmp(sel,'FULL')
    out=3*f1*trQ;
end
end

function [f1 f2]=compute_eshelby_factors(E,nu_c,nu_t)
%could in principle compute non-shperical or whatever, but sticking with
%this for now

% The following two lines are correct for using biphasic greens function
% but for the purposes of this plot we use the instantaneous response
%f1=((1-nu_t)/(1+nu_t))./(1+2.*E*(1-2*nu_c)/(1+nu_t));
%f2=(15/4*(1-nu_t)/(4-5*nu_t))./(1+E/2*(1+nu_c)*(7-5*nu_t)/(1+nu_t)/(4-5*nu_t));

f1=1./(3*(1+2.*E*(1-2*nu_c)/(1+nu_t)));
f2=(1/2)./(1+E/2*(1+nu_c)*(7-5*nu_t)/(1+nu_t)/(4-5*nu_t));
end


function ret=eshelby_mismatch_in(Efrac,nu_c,nu_t)
%returns trace value since that's all we need
if nargin<3
    nu_t=nu_c;
end
if nargin<2
    nu_t=0.4;nu_c=0.4;
end
E=Efrac; %this should be E_tissue/E_cell
ret=(1+(E.*(1-2*nu_c)/(1-2*nu_t)-1)./(2*(1-nu_c)/(1+nu_t).*E + 1));
end