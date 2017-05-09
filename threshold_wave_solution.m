function [vel veld G]=threshold_wave_solution(fig,G)

global N;
global E;
global E0;
global r;
global k;
global a;
global b;
global delta;
global modtype;
global stype;
global htype;

if nargin<1
    fig=1;
end
if nargin<2
    G=linspace(0.1,5,200);
end
r=50;N=10;
k=1.5;a=0.9;b=0.10;delta=1.4;
modtype=0;stype=0;htype=0; %d=1; no strain modulation this "less realistic" stuff

%r=0.01;gamma=3;N=6;
%k=8;a=0.75;b=0.25;delta=1.0;

vel=zeros(length(G),1);veld=zeros(length(G),1);
strain=zeros(length(G),1);straind=zeros(length(G),1);
%E0=37.5; %E0=75; %Larry Taber's microindentation measurements
%E0=200; %E0=100; %Stephanie's micropipette aspiration measurements
%E0=55;
E0=55;t0=0.001*[1 1]';%for testing purposes
for i=1:length(G)
    E=E0*G(i);
    [t td]=wave_solution(t0);
    %[t td]=wave_solution(t0,fig+2);
    if(~isnan(t))
        vel(i)=1./t;
        t0(1)=t*0.9;
        veld(i)=1./td;
        strain(i)=F(E,modtype)/90.*vel(i)./(2*E);
    end
    if(~isnan(td))
        veld(i)=1./td;
        t0(2)=td*0.9;
        straind(i)=F(E,modtype)/90.*veld(i)./(2*E);
    end
end
[es ese]=expt_strain_vals(1);
[ev eve]=expt_vel_vals(1);
figure(fig);close(fig);figure(fig);subplot(2,1,1);%plot(G,vel,'r.');
hold on;plot(G,vel/100,'r.');errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'bs');
xlabel('E/E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
subplot(2,1,2);hold on;plot(G,strain,'r.');errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'bs');
xlabel('E/E^*');ylabel('\epsilon (% strain)');hold off;
subplot(2,1,1);hold off;
figure(fig+1);close(fig+1);figure(fig+1);subplot(2,1,1);%plot(G,vel,'r.');
hold on;plot(G,veld/100,'r.');errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'bs');
xlabel('E/E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
subplot(2,1,2);hold on;plot(G,straind,'r.');errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'bs');
xlabel('E/E^*');ylabel('\epsilon (% strain)');hold off;
subplot(2,1,1);hold off;
figure(100);semilogx(G,F(E0*G,modtype),'r-');
ofp=fopen('velnum.dat','w');
A=[G',vel/100,strain];
fprintf(ofp,[repmat('%f ', 1, size(A, 2)), '\b\n'], A.');
fclose(ofp);

function [tval tvald]=wave_solution(t0,fig)

dt=logspace(-4,0,300);
if nargin<1
    t0=[0.005 0.005];
end

%monopole excitation
left=lhs(dt);
right=rhs(dt);
tval=fzero(@fun,t0(1));

%dipole excitation
leftd=lhs_dipole(dt);
rightd=rhs_dipole(dt);
tvald=fzero(@fun_dipole,t0(2));

if nargin>1
%     figure(fig);
%     axis([min(dt),max(dt),-left(1),2*left(1)])
%     plot(log(dt),log(left),'r--');hold on;
%     plot(log(dt),log(right),'b--');
%     plot(log(dt),log(left-right),'g--');
%     plot(log(tval),log(rhs(tval)),'ks');
%     hold off;
    figure(fig+1);
    %axis([min(dt),max(dt),-leftd(1),2*leftd(1)])
    plot(log(dt),log(leftd),'r--');hold on;
    plot(log(dt),log(rightd),'b--');
    plot(log(dt),log(leftd-rightd),'g--');
    plot(log(tvald),log(rhs_dipole(tvald)),'ks');
    hold off;
end

function result=fun(x)
result=rhs(x)-lhs(x);

function result=lhs(~)
global E;global r;global modtype;
result=r.*sqrt(16*pi*E)./(F(E,modtype).^2);

function result=rhs(x)
global E;global N;
if N==0
    result=x.^(-3/2).*polylog(0.5,exp(-1./(4.*E.*x)));
else
    result=0;
    for n=1:N
        result=result+sqrt(1./(n.*(x.^3))).*(exp(-n./(4.*E.*x)));
    end
end

function result=fun_dipole(x)
result=rhs_dipole(x)-lhs_dipole(x);

function res=rhs_dipole(x)
global E;global N;
if N==0
    %polylog solution
else
    res=0;
    for n=1:N
        res=res+x.^(-3/2).*exp(-n./(4.*E.*x)).*(1./(2.*E.*x.*sqrt(n))-1./(n.^(3/2)));
    end
end

function res=lhs_dipole(~)
global E;global r;global modtype;
res=r.*sqrt(16*pi*E)./(F(E,modtype).^2);

function mod=F(E,type)

global k;
global a;
global b;
global delta;
global E0;
global stype;

e=E/E0;

extramod=softfcn(e,stype);

if nargin<2
    type=1;
end

e=e-delta; %shift
%mod=1;
if type==1
    if e>0
        mod=(a*exp(-2*k.*e)+b)./(exp(-2*k.*e)+1);
    else
        mod=(a+b*exp(2*k.*e))./(1+exp(2*k.*e));
    end
elseif type==2
    if e>0
        mod=b*ones(size(E));
    else
        mod=a*ones(size(E));
    end
elseif type==0
    mod=ones(size(E)); %no modification
end
mod=mod.*extramod;

function mod=softfcn(e,type)

I=ones(size(e));
mod=I;
if type>0
    e=e-0.3.*I;
    c=1.5;
    mod=I./(exp(-2.*c.*e)+I);
end

function mod=H(E,type)
global E0;

e=E/E0;

if type==0
    mod=1; %no modification
    return;
elseif type==1
    if e<1
        mod=e;
    else
        mod=1;
    end
end

function [val err]=expt_strain_vals(bool)
if bool==0
    val=0;err=0;
else
    val=[0.2 0.5 0.7 1 1.5 2.5;0.125*[0.27 0.2 0.72 1.0 0.8 0.5]];
    err=.125*[.15 .12 .22 .03 .04 .04;.14 .13 .25 .06 .05 .05];
end

function [val err]=expt_vel_vals(bool)
if bool==0
    val=0;err=0;
else
    val=[0.2 0.5 0.7 1 1.5 2.5;2 10.5 13.5 21.5 26 0];
    err=[2 2 0 5.5 5 0;3 4 1.5 6 5.5 0];
end