function [Ef,nuf] = compositeTissueCalc(Elist,fig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
    Elist=linspace(0.1,3,30);
end
if nargin<2
    fig=[];
end

Ec=0.47;
Ef=zeros(size(Elist));nuf=zeros(size(Elist));
%[K,G]=calculateKGfromENu(1,0.4); %kid tested, mother approved
%[E,nu]=calculateENufromKG(K,G); %kid tested, mother approved
for i=1:numel(Elist)
    Evec=[Elist(i);Ec];nuvec=[0.4;0.4];Vvec=[1-4*pi/24;4*pi/24];
    [Kvec,Gvec]=calculateKGfromENu(Evec,nuvec);
    [Kf,Gf]=compositeEshelbyMaterial(Kvec,Gvec,Vvec);
    [Ef(i),nuf(i)]=calculateENufromKG(Kf,Gf);
end
if ~isempty(fig)
    figure(fig);plot(Elist,Ef,'r-');
end
end

function [K,G]=calculateKGfromENu(E,nu)
K=E./(3*(1-2*nu));G=E./(2*(1+nu));
end

function [E,nu]=calculateENufromKG(K,G)
E=9*K.*G./(3*K+G);nu=(3*K-2*G)./(2*(3*K+G));
end

function [Kf,Gf]=compositeEshelbyMaterial(Kvec,Gvec,Vvec)
if ~isequal(size(Kvec),size(Gvec)) || ~isequal(size(Gvec),size(Vvec)) || ~isequal(size(Kvec),size(Vvec))
    error('');
end
K0=Kvec(1);G0=Gvec(1);
Ki=Kvec(2:end);Gi=Gvec(2:end);Vi=Vvec(2:end);
[~,nui]=calculateENufromKG(Ki,Gi);

Gf=G0/(1+sum(Vi.*(1-Gi/G0))/(1+2*sum((4-5*nui)./(15-15*nui).*(Gi/G0-1))));
Kf=K0/(1+sum(Vi.*(1-Ki/K0))/(1+1/3*sum((1+nui)./(1-nui).*(Ki/K0-1))));
end