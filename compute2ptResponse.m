function [ output_args ] = compute2ptResponse(fig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

linplotflag=1;
logplotflag=1;

if nargin<1
    fig=1;
end
%G=linspace(0,100,1000);
G=logspace(-0.8,0.5,100);
Tr_eps=zeros(size(G));
nu=0.4;Q0=0.1;Ec=0.47;
cellQ=[2 0 0;0 -2*nu 0;0 0 -2*nu];

r=[1;0;0];r2=r'*r;

Tin=eshelby_mismatch_in(G/Ec,nu,nu);G1=zeros(size(G));G2=zeros(size(G));
[f1,f2]=compute_eshelby_factors(G/Ec,nu,nu);

for i=1:numel(G)
    E=G(i);F1=f1(i);F2=f2(i);
    D=E*(1-nu)/((1+nu)*(1-2*nu));
    epstemp=eshelby_mismatch_out(E/Ec,Q0*cellQ,nu,nu,F1,F2);
    %epstemp=eshelby_mismatch_out(E/Ec,Q0*cellQ,nu,nu);
    Q=E/(1+nu)*(epstemp+nu/(1-2*nu)*eye(3)*trace(epstemp)); %convert to matrix stress
    trQ=trace(Q);projQ=r'*(Q*r)/r2;
    G1(i)=G1Static(r,D,trQ);G2(i)=G2Static(r,D,projQ);
    Tr_eps(i)=-Tin(i)*(G1(i)-G2(i));
end

%import Jason's finite elt data
%fEltData=csvread('./finiteEltComsolData/cell_interaction_static.csv',1);
scriptCalculateCompositeModulusFactor; %use a script which contains newer data
%now these variables exist

if linplotflag==1
    figure(fig);
    %Emod=compositeTissueCalc(fEltData(:,1));%Emod=fEltData(:,1); %base case
    plot(1.6*G,Tr_eps,'b-');hold on;%plot(fEltData(:,1),2*fEltData(:,2),'b--',Emod,2*fEltData(:,2),'b-');hold off;legend('Analytics','Finite Elt','Composite finite elt','Location','NorthEast');
    plot(Em,str,'r.');hold off;legend('Point CM','Finite CM','Location','NorthEast');
    axis([0 5 0 0.02]);xlabel('Tissue Young''s modulus (kPa)');ylabel('Strain in neighboring cell');
end
if logplotflag==1
    figure(fig+1);
    %Emod=compositeTissueCalc(fEltData(:,1));%Emod=fEltData(:,1); %base case
    loglog(1.6*G,Tr_eps,'b-');hold on;%loglog(fEltData(:,1),2*fEltData(:,2),'b--',Emod,2*fEltData(:,2),'b-');hold off;legend('Analytics','Finite Elt','Composite finite elt','Location','NorthEast');
    loglog(Em,str,'r.');hold off;legend('Point CM','Finite CM','Location','NorthEast');
    axis([0 5 0 0.02]);xlabel('Tissue Young''s modulus (kPa)');ylabel('Strain in neighboring cell');
end

%figure(fig+1);plot(G,Tin,'r.');
%figure(fig+2);plot(G,f1,'r.');
%figure(fig+3);plot(G,Tin.*f1,'r.');
%figure(fig+4);plot(G,G1,'r.',G,G2,'b.');
plot2svg('./ComparePointFiniteCMExcite_linplot.svg',1);
plot2svg('./ComparePointFiniteCMExcite_logplot.svg',2);
end

function out=G1Static(r,D,trQ)
r3=(r'*r)^(3/2);
out=trQ./(4*pi*D*r3);
end

function out=G2Static(r,D,projQ)
r3=(r'*r)^(3/2);
out=3.*projQ./(4*pi*D*r3);
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
if nargin<4
    nu_t=nu_c;
end
if nargin<2
    nu_t=0.4;nu_c=0.4;
end
if nargin<5
    %only compute these if they don't already exist
    [f1 f2]=compute_eshelby_factors(E,nu_c,nu_t);
end
ret=f1*eye(3)*TrQ+f2*(Q+Q'-2/3*eye(3)*TrQ);
end

function [f1 f2]=compute_eshelby_factors(E,nu_c,nu_t)
%could in principle compute for non-shperical inclusions, but sticking with
%this for now
f1=((1-nu_t)/(1+nu_t))./(1+2.*E*(1-2*nu_c)/(1+nu_t));
f2=(15/4*(1-nu_t)/(4-5*nu_t))./(1+E/2*(1+nu_c)*(7-5*nu_t)/(1+nu_t)/(4-5*nu_t));
end