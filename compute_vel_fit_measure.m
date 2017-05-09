function lsqrs=compute_vel_fit_measure(E_v,dataopt,fitopt)
%computes a velocity fit
if nargin<2
    dataopt=0;
end
if nargin<3
    fitopt=0;
end
if dataopt==1
    val=[0.2 0.5 0.7 1 1.5;2 10.5 13.5 21.5 26];
elseif dataopt==0
    val=[0.2 0.5 0.7 1 1.5 2.5;2 10.5 13.5 21.5 26 0];
    %err=[2 2 0 5.5 5 0;3 4 1.5 6 5.5 0];
elseif dataopt==2
    val=[0.5 0.7 1 1.5;2 10.5 13.5 21.5];
else
    error('bad option!');
end
%compute least squares measure
E=E_v(:,1);v=E_v(:,2);
if fitopt==1
    lsqrs=least_squares(val,[E,v]);
elseif fitopt==0
    lsqrs=least_squares_xdiv(val,[E,v],20);
end
end

function lsq=least_squares(expt_EV,sim_EV)
sze=size(expt_EV);szs=size(sim_EV);
if sze(1)==2
    expt_EV=expt_EV';
elseif sze(2)~=2
    error('Oops!\n');
end
if szs(1)==2
    sim_EV=sim_EV';
elseif szs(2)~=2
    error('Oops!\n');
end
v1=expt_EV(:,2);E1=expt_EV(:,1);
v2=sim_EV(:,2);E2=sim_EV(:,1);
idx=zeros(length(v1),1);
for i=1:length(E1)
    [~,idx(i)]=min(abs(E2-E1(i))); %locate closest E value to i'th entry
end
vec=v1-v2(idx); %use linear least squares at closest value
lsq=vec*vec';
end

function res=least_squares_xdiv(expt_EV,sim_EV,ratio)
%here v1 and v2 have both E and v components
sze=size(expt_EV);szs=size(sim_EV);
if sze(1)==2
    expt_EV=expt_EV';
elseif sze(2)~=2
    error('Oops!\n');
end
if szs(1)==2
    sim_EV=sim_EV';
elseif szs(2)~=2
    error('Oops!\n');
end
if nargin<3
    ratio=20;
end
% if ratio<10
%     error('ratio < 10 : odd least squares metric');
% end
v1=expt_EV(:,2);E1=expt_EV(:,1);
v2=sim_EV(:,2);E2=sim_EV(:,1);
lsq=zeros(length(v1),1);
for i=1:length(v1)
    %minimum distance using an ellipsoidal metric
    lsq(i)=min((v1(i)-v2).^2+ratio^2*(E1(i)-E2).^2);
    if i==length(v1)
        lsq(i)=lsq(i)*1;
    end
end
res=sum(lsq);
end

