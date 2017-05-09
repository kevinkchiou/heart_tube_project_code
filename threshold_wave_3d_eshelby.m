function [v strain converge] = threshold_wave_3d_eshelby(fig,G,alphabar,drag,mod,eshelby_params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global N;
global E0;
global E;
global r; %alpha / Q
global nu;
global Gamma;
global modtype;
global Qtype;
global conv;

%eshelby factors to speed up computation
global F1;
global F2;

%fits to cell-on-gel data
global Esat;
global Ec;
global eps_trans;

if nargin<1
    fig=1;
end
if nargin<2
    G=linspace(0.1,3); %include expt vals
end
%E0 = 1.3kPa for E4 embryos.
if nargin<3
    %r=1.47;Gamma=0.0015;Q=0.10;modtype=1;Esat=1;Ec=1;eps_trans=1;%constant stress source (eshelby)
    %r=1.40;Gamma=0.0012;Q=0.10;modtype=1;Esat=0.7323;Ec=1.2727;eps_trans=0.904;%constant stress source (eshelby)
    %r=0.2;Gamma=0.006;Q=0.10;modtype=2;Esat=0.7323;Ec=1.2727;eps_trans=0.904;%constant stress source (eshelby)
    N=5;r=1.25;Gamma=0.015;Q=0.075;modtype=2;Esat=0.7323;Ec=0.7;eps_trans=1.0;%constant stress source (eshelby)
    %r=0.4;Gamma=0.006;Q=0.1;modtype=1;
    %r=0.4;Gamma=0.005;Q=0.2;modtype=2;%stall force stress source (eshelby)
    %r=0.2;Gamma=0.008;Q=0.1;modtype=2;%stall force stress source (eshelby) testing
    %N=5;r=0.6;Gamma=0.0017;Q=0.075;modtype=2;Esat=0.7323;Ec=0.7;eps_trans=1.0;%constant stress source (eshelby)
    r=0.5377;Gamma=0.0029;Q=0.10;Esat=0.6823;Ec=0.5947;eps_trans=1;
    r=0.65; %near threshold limit of steady state solution given parameters right above.
    %N=20;r=4.5;Gamma=0.0008;Q=0.15;modtype=2;Esat=0.7323;Ec=0.7;eps_trans=1.05; %testing
else
    r=alphabar;
end
if nargin==3
    Gamma=0.011; %guessing
end
if nargin==4
    Gamma=drag;
end
if nargin==5
    modtype=mod;
end
if nargin==6
    r=alphabar;
    Gamma=drag;
    modtype=mod;
    Q=eshelby_params(1);
    Esat=eshelby_params(2);
    Ec=eshelby_params(3);
    eps_trans=eshelby_params(4);
end

%parameters fit to Stephanie's data given a couple of fits. A stress
%saturation model (Esat), Cell stiffness of Eshelby inclusion (Ec), and
%an overall scale of strain transmission (eps_trans). These were fit with
%the function cell_stiffness_fit_residual(). Can do this for Engler's data
%as well

%Esat=0.7323;Ec=1.2727;eps_trans=0.9040;


G=sort(G);
%[G,tidx]=sort(G); %include expt values and sort
%invidx(tidx)=1:length(G); %inverse index to compare to expt velocities

%N=5;
N=10;nu=0.4;Qtype=1;
%modtype=2; %1 is constant stress source, 2 is linear stress until stall,

if Qtype==1
    Q0=[2 0 0;0 -2*nu 0;0 0 -2*nu];
else
    Q0=[2 0 0;0 -1 0;0 0 -1];
end
vel=zeros(length(G),1);strain=zeros(length(G),1);converge=zeros(length(G),1);
%E0=37.5; %E0=75; %Larry Taber's microindentation measurements
%E0=200; %E0=100; %Stephanie's micropipette aspiration measurements
E0=1;t0=0.005*[1 1]'; %t0=0.005*[1 1]'; %for delta(t) solutions
[f1 f2]=compute_eshelby_factors(E0*G/Ec,nu,nu);
for i=1:length(G)
    %global variables to make computation quicker
    E=E0*G(i);F1=f1(i);F2=f2(i);
    t=wave_solution_theta(t0);
    converge(i)=conv; %convergence test using global variable 
    %t=wave_solution_theta(t0,20); %use to diagnose root finding
    if(~isnan(t))
        vel(i)=1./t;
        t0(1)=t*0.9;
        tempQ=E*eps_trans*eshelby_mismatch_out(E/Ec,Q0,nu,nu,F1,F2);
        %tempQ=E*eps_trans*eshelby_mismatch_out(E/Ec,[1 0 0;0 -nu 0;0 0 -nu],nu,nu,F1,F2);
        %baseline physical strain
        D=E*(1-nu)/(Gamma*(1+nu)*(1-2*nu));
        strain(i)=Q*tempQ(1,1)/(D^2*Gamma).*vel(i);
    end
end
if modtype==2
    fmod_linear=min([G/Esat;ones(1,numel(G))]);
    strain=strain.*fmod_linear';
elseif modtype==3
end
[es ese]=expt_strain_vals(2);
[ev eve]=expt_vel_vals(2);
[ec ece]=expt_cellgel_vals(1);
if fig>0
    figure(fig);close(fig);set(figure(fig),'Position',[0 0 600 1000]);subplot(2,1,1);%plot(G,vel,'r.');
    hold on;plot(G,vel/50,'r.');
    errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'bs');
    xlabel('E / E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
    subplot(2,1,2);hold on;plot(G,strain,'r.');
    errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'bs');
    xlabel('E / E^*');ylabel('\epsilon (% strain)');hold off;
    subplot(2,1,1);hold off;
end
if (modtype==1 && fig>0)
    Gother=linspace(0.1,45,200);
    figure(100);
    L1=eps_trans*compute_eshelby_factors(E0*Gother/Ec,nu,nu);
    %L2=eshelby_mismatch_in(E0*Gother,nu,nu);
    plot(log10(Gother),3*L1,'r-');hold on; %factor of three when trace is taken
    %set(plot(log10(Gother),L2,'k-'),'LineWidth',2); %no factor of three - trace is included
    errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
    %title('Constant stress source, spherical Eshelby cell');
    xlabel('log_{10} ( E / E^* )');
    %legend('Eshelby out','Eshelby in','Experiment','Location','NorthWest');
    %legend('Eshelby out','Experiment','Location','NorthEast');
    hold off;
elseif (modtype==2 && fig>0)
    Gother=linspace(0.1,45,200);
    figure(101);
    L1=eps_trans*compute_eshelby_factors(E0*Gother/Ec,nu,nu);
    L2=eshelby_mismatch_in(E0*Gother/Ec,nu,nu);
    %following line is linear dependence of cell stress activation on E
    fmod_linear_Gother=min([Gother/Esat;ones(1,numel(Gother))]);
    plot(log10(Gother),3*L1.*fmod_linear_Gother,'r-');hold on; %factor of three when trace is taken
    set(plot(log10(Gother),L2,'k-'),'LineWidth',2);hold on; %no factor of three - trace is included
    errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
    %title('Linear stress until stall, spherical Eshelby cell');
    xlabel('log_{10} ( E / E^* )');
    %legend('Linear stress until stall','Eshelby in','Experiment','Location','NorthWest');
    hold off;
elseif modtype==3 && fig>0
    Gother=linspace(0.1,45,200);
    figure(101);
    plot(log10(Gother),M(E0*Gother,modtype)/.75,'r-');hold on;
    set(plot(log10(G),M(E0*G,modtype)/0.75,'k-'),'LineWidth',2);
    errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
    hold off;
end
%ofp=fopen('velnum.dat','w');
%A=[G',vel/50,strain];
%fprintf(ofp,[repmat('%f ', 1, size(A, 2)), '\b\n'], A.');
%fclose(ofp);
v=vel/50; %output velocity in physical units
end

%theta is a H(t)H(\tau - t) activation instead of delta(t)
function tval=wave_solution_theta(t0,fig)

dt=logspace(-4,0,300); %row vector

%3d dipole excitation
left=ones(1,length(dt)); %row vector
right=rhs_theta(dt);
if nargin<1
    t0=0.01;
end
t03d=t0(1)*.85;
opts = optimset('Display','off');
%tval=fzero(@fun_theta,[1e-6 0.1],opts);
tval=fzero(@fun_theta,t03d,opts);

if nargin>1
    figure(fig);
    plot(log10(dt),log10(left),'r--');hold on;
    plot(log10(dt),log10(right),'b--');
    plot(log10(dt),log10(left-right),'g--');
    plot(log10(tval),log10(rhs_theta(tval)),'ks');
    axis(log10(2.71)*[-10 -4 -8 8]);
    hold off;
end
end

function result=fun_theta(x)
result=rhs_theta(x)-1;
end

function res=W1(a,a2,trQ,idx)
if trQ==0
    res=zeros(size(a2));return;
end
if sum(a2<0)>0
    res(a2<0)=NaN;
end
if nargin<4
    idx=ones(size(a));
end
res(idx==0)=0;bm=(idx==1 & a2>0);
res(bm)=trQ*(erfc(a(bm))+2*a(bm)/sqrt(pi).*exp(-a2(bm)));
end

function res=W2(a,a2,projQ,idx)
if projQ==0
    zeros(size(a2));return;
end
if sum(a2<0)>0
    res(a2<0)=NaN;
end
if nargin<4
    idx=ones(size(a));
end
res(idx==0)=0;bm=(idx==1 & a2>0);
res(bm)=projQ*(3*erfc(a(bm))+(0.5*(2*a(bm)).^3 + 6*a(bm))/sqrt(pi).*exp(-a2(bm)));
end

function result=rhs_theta(x)
global E;global N;global Gamma;global nu;global r;global F1;global F2;global E0;
global modtype;
global Ec;global Esat;global eps_trans;
global Qtype;
global conv; %convergence test
if sum(isnan(x))>0
    result=zeros(size(x));
    result(isnan(x))=NaN;return;
end

if  Qtype==1
    Q0=[2 0 0;0 -2*nu 0; 0 0 -2*nu];
else
    Q0=[2 0 0;0 -1 0;0 0 -1];
end

D = E.*(1-nu)/(Gamma*(1+nu)*(1-2*nu));
y = 1./(4*D*x); %hopefully will fix complaints
coeff = (1+nu)/(4*pi*(1-nu)*r); %simplified from E/(Gamma(1-2nu)*D) - stress trace GF coeff.
res=0;tau=1e-1; %contraction time
if modtype==1
    cellQ=Q0; %constant stress source
elseif modtype==2
    cellQ=min(E/E0/Esat,1)*Q0; %linear growth in excitation
end
%Take cellQ as cell strain
epstemp=eps_trans*eshelby_mismatch_out(E/E0/Ec,cellQ,nu,nu,F1,F2); %transmit to matrix strain
Q=E*(1/(1+nu)*epstemp+nu/(1-2*nu)*eye(3)*trace(epstemp)); %convert to matrix stress
trQ=trace(Q);
%now compute Green function terms for myocyte array - this is matrix stress
%propagation using G_ijkl across various distances.
for nx=1:N+1
    %trQ=0;projQ=2; %old excitation
    dr=[nx;0;0];dr2=nx^2;projQ=dr'*(Q*dr)/dr2;%new excitation in purely x direction
    a2=nx.*y;a=sqrt(a2);
    res=res+coeff.*(W1(a,a2,trQ)+W2(a,a2,projQ))/(dr2)^(3/2);
    %need to add component that turns them off
    offidx=nx.*x - tau > 0; %which myocytes are turned off
    a2_d=nx.^2./(4*D*(nx.*x - tau));a_d=sqrt(a2_d);
    res=res-coeff.*(W1(a_d,a2_d,trQ,offidx)+W2(a_d,a2_d,projQ,offidx))/(dr2)^(3/2);
    for ny=1:N+1 %symmetric under exchange ny <-> -ny gives factor of 2
        a2 = (nx+(ny^2/nx)).*y;a=sqrt(a2);
        %trQ=0;projQ=(2*nx.^2-ny.^2)/(nx.^2+ny.^2); %old excitation
        dr=[nx;ny;0];dr2=nx^2+ny^2;projQ=(dr'*Q*dr)/dr2; %new excitation
        res=res+coeff.*(W1(a,a2,trQ)+W2(a,a2,projQ))/(dr2)^(3/2);
        %add off component via manual heaviside
        a2_d=(nx.^2+ny.^2)./(4*D*(nx.*x - tau));a_d=sqrt(a2_d);
        res=res-coeff.*(W1(a_d,a2_d,trQ,offidx)+W2(a_d,a2_d,projQ,offidx))/(dr2)^(3/2);
    end
    if nx==N
        NMinusOneVal=res;
    end
    if nx==(N+1)
        NVal=res;
    end
end
convergenceTest=(NVal-NMinusOneVal)/NVal;conv=convergenceTest;%test for convergence

result_temp=res*(1-2*nu)/E; %convert from Tr(s_m) into Tr(eps_m)
fac=eshelby_mismatch_in(E/E0/Ec,nu,nu); %eshelby stiffness mismatch - induced
result=result_temp.*fac; %transmit Tr(eps_m) into Tr(eps_c)
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

function [ret f1 f2]=eshelby_mismatch_out(Efrac,Q,nu_c,nu_t,f1,f2)
%this can't be just a trace - this must be a full tensor which gives us our
%Q^*_ij in the tissue after taking in a Q_kl from the cell.
E=Efrac;TrQ=trace(Q); %E_tissue/E_cell, and trace of the cell excitation
if nargin<3
    nu_t=nu_c;
end
if nargin<2
    nu_t=0.4;nu_c=0.4;
end
if nargin<4
    %only compute these if they don't already exist
    [f1 f2]=compute_eshelby_factors(E,nu_c,nu_t);
end
ret=f1*eye(3)*TrQ+f2*(Q+Q'-2/3*eye(3)*TrQ);
end

function [f1 f2]=compute_eshelby_factors(E,nu_c,nu_t)
%could in principle compute non-shperical or whatever, but sticking with
%this for now
f1=((1-nu_t)/(1+nu_t))./(1+2.*E*(1-2*nu_c)/(1+nu_t));
f2=(15/4*(1-nu_t)/(4-5*nu_t))./(1+E/2*(1+nu_c)*(7-5*nu_t)/(1+nu_t)/(4-5*nu_t));
end
