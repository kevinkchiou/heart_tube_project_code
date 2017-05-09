function [vel vel_strain strain straind lsqval]=threshold_wave_3d_solution(fig,exptopt,G,alphabar,drag)

global N;
global E;
global E0;
global r; %alpha / Q
global nu;
global Gamma;
global k;
global a;
global b;
global delta;
global modtype;

if nargin<1
    fig=1;
end
if nargin<2
    exptopt=0;
end
if nargin<3
    G=linspace(0.1,3); %include expt vals
end
Q=1/2200; %just as a start adjust as necessary
if nargin<4
    %r=70;Gamma=0.05;Q=1/2200; %velocity fit with T^2
    %r=1350.9;Gamma=0.0344;Q=1/2200 %upper cutoff fit with T^2
    %r=500;Gamma=0.05; %upper cutoff guess
    %r=256;Gamma=0.0331%strain activation F^2 fit with upper cutoff
    %r=.01;Gamma=0.0963%F^1 stress activation fit without upper cutoff
    %r=268.3;Gamma=0.0331;%F^1 strain activation fit with upper cutoff
    %r=281;Gamma=0.0331;%F^2 strain activation fit without upper cutoff
    %r=0.121;Gamma=0.087;Q=1/3700;%F^1 strain activation fit without upper cutoff
    
    %r=0.23;Gamma=0.01;Q=1/530; %heaviside (good) guess T^2 with upper cutoff
    r=0.1326;Gamma=0.0165;Q=1/630; %heaviside fit with T^1 no upper cutoff
else
    r=alphabar;
end
if nargin==4
    %Gamma=0.05; %velocity fit
    %Gamma=0.0344; %upper cutoff fit from parameter_best_fit_search()
    %Gamma=0.05; %upper cutoff fit
    %Gamma=0.0331;%strain activation F^2 fit with upper cutoff
    %Gamma=0.0963;%F^1 stress activation fit without upper cutoff
    %Gamma=0.0331;%F^1 strain activation fit with upper cutoff
    %Gamma=0.0331;%F^2 strain activation fit without upper cutoff
    %Gamma=0.087;%F^1 strain activation fit without upper cutoff
    
    Gamma=0.011; %guessing
end
if nargin==5
    Gamma=drag;
end
if exptopt==1
    G=[G, 0.2, 0.5, 0.7, 1, 1.5 2.5];
end
[G,tidx]=sort(G); %include expt values and sort
invidx(tidx)=1:length(G); %inverse index to compare to expt velocities
N=5;nu=0.4;
k=1.5;a=0.9;b=0.10;delta=1.4;
modtype=2;
%k=8;a=0.75;b=0.25;delta=1.0;

t1suppress=0; %suppresses plots of d=1 cases
vel=zeros(length(G),1);vel_strain=zeros(length(G),1);
strain=zeros(length(G),1);straind=zeros(length(G),1);
%E0=37.5; %E0=75; %Larry Taber's microindentation measurements
%E0=200; %E0=100; %Stephanie's micropipette aspiration measurements
E0=1;t0=0.005*[1 1]'; %t0=0.005*[1 1]'; %for delta(t) solutions
for i=1:length(G)
    E=E0*G(i);
    %[t td]=wave_solution(fig+2);
    %[t t_strain]=wave_solution(t0);t1suppress=1;
    %[t t_strain]=wave_solution(t0,10); %use to diagnose root finding
    [t t_strain]=wave_solution_theta(t0);t1suppress=1;
    %[t t_strain]=wave_solution_theta(t0,20);t1suppress=0; %use to diagnose root finding
    if(~isnan(t))
        vel(i)=1./t;
        t0(1)=t*0.9;
        %km=vel(i)./(2*E).*(-1+sqrt(1+gamma*E./vel(i)^2));
        %kp=vel(i)./(2*E).*(-1-sqrt(1+gamma*E./vel(i)^2));
        if modtype==1
            strain(i)=L(E)/2000.*vel(i).^(7/2)./(2*(1-2*nu)*E/Gamma);
            strain(i)=L(E)/2000.*vel(i)./(2*(1-2*nu)*E/Gamma);
        elseif modtype==2
            %strain(i)=F(E)*Q.*vel(i).^(7/2)./(2*(1-2*nu)*E/Gamma);
            strain(i)=F(E)*Q.*vel(i)./(2*(1-2*nu)*E/Gamma);
        elseif modtype==3
            %strain(i)=M(E)/1000.*vel(i).^(7/2)./(2*(1-2*nu)*E/Gamma);
            strain(i)=M(E)/1000.*vel(i)./(2*(1-2*nu)*E/Gamma);
        end
%         if(i==floor(length(G)/2))
%             figure(101);xp=linspace(-5*floor(1/km),0);xn=linspace(0,3);
%             plot([xp xn],[exp(km*xp) exp(kp*xn)],'b-');
%         end
    end
    if(~isnan(t_strain))
        vel_strain(i)=1./t_strain;
        t0(2)=t_strain*0.9;
        if modtype==1
            %straind(i)=L(E)/120.*vel_1d(i)./(2*E/Gamma);
            straind(i)=L(E)/2000.*vel(i)./(2*(1-2*nu)*E/Gamma);
        elseif modtype==2
            %straind(i)=F(E)/90.*vel_1d(i)./(2*E/Gamma);
            straind(i)=F(E)/1200.*vel(i)./(2*(1-2*nu)*E/Gamma);
        elseif modtype==3
            %straind(i)=M(E)/90.*vel_1d(i)./(2*E/Gamma);
            straind(i)=M(E)/1000.*vel(i)./(2*(1-2*nu)*E/Gamma);
        end
    end
end
[es ese]=expt_strain_vals(1);
[ev eve]=expt_vel_vals(1);
[ec ece]=expt_cellgel_vals(1);
if fig>0
    figure(fig);close(fig);figure(fig);subplot(2,1,1);%plot(G,vel,'r.');
    hold on;plot(G,vel/100,'r.');errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'bs');
    xlabel('E/E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
    subplot(2,1,2);hold on;plot(G,strain,'r.');errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'bs');
    xlabel('E/E^*');ylabel('\epsilon (% strain)');hold off;
    subplot(2,1,1);hold off;
end
if t1suppress==0 && fig>0
    figure(fig+1);close(fig+1);figure(fig+1);subplot(2,1,1);%plot(G,vel,'r.');
    hold on;plot(G,vel_strain/100,'r.');errorbar(ev(1,:),ev(2,:),eve(1,:),eve(2,:),'bs');
    xlabel('E/E^*');ylabel('v (mm/sec)');hold off;%axis([0 max(G) 0 30])
    subplot(2,1,2);hold on;plot(G,straind,'r.');errorbar(es(1,:),es(2,:),ese(1,:),ese(2,:),'bs');
    xlabel('E/E^*');ylabel('\epsilon (% strain)');hold off;
    subplot(2,1,1);hold off;
end
if modtype==1 && fig>0
    Gother=linspace(0.1,30,200);
    figure(100);
    plot(log10(Gother),L(E0*Gother),'r-');hold on;
    set(plot(log10(G),L(E0*G),'k-'),'LineWidth',2);
    errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
    hold off;
elseif modtype==2 && fig>0
    Gother=linspace(0.1,30,200);
    figure(101);
    plot(log10(Gother),F(E0*Gother),'r-');hold on;
    set(plot(log10(G),F(E0*G),'k-'),'LineWidth',2);
    errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
    hold off;
elseif modtype==3 && fig>0
    Gother=linspace(0.1,30,200);
    figure(101);
    plot(log10(Gother),M(E0*Gother,modtype)/.75,'r-');hold on;
    set(plot(log10(G),M(E0*G,modtype)/0.75,'k-'),'LineWidth',2);
    errorbar(log10(ec(1,:)),ec(2,:),ece(1,:),ece(2,:),'bs');
    hold off;
end
ofp=fopen('velnum.dat','w');
A=[G',vel/100,strain];
fprintf(ofp,[repmat('%f ', 1, size(A, 2)), '\b\n'], A.');
fclose(ofp);
vel=vel/100; %output velocity in physical units
vel_strain=vel_strain/100;lsqval=NaN;
if exptopt==1
    %v=vel(invidx(end-5:end));lsqval=compute_vel_fit_measure([G(invidx(end-5:end))',v],0);
    v=vel(invidx(end-5:end-1));lsqval=compute_vel_fit_measure([G(invidx(end-5:end-1))',v],1);
    %v=vel(invidx(end-4:end-1));lsqval=compute_vel_fit_measure([G(invidx(end-4:end-1))',v],2);
end

function [tval tval_strain]=wave_solution(t0,fig)

dt=logspace(-5,0,300); %row vector

%3d dipole excitation
left=ones(1,length(dt)); %row vector
right=rhs(dt);
if nargin<1
    t0=[0.01 0.01];
end
if length(t0)<2
    t03d=t0;t01d=t0;
    t0_strain=t0;
else
    t03d=t0(1)*.75;t01d=t0(2)*.75;
    t0_strain=t0(2)*.75;
end
left_strain=ones(1,length(dt));
right_strain=rhs_strain(dt);
opts = optimset('Display','off');
tval=fzero(@fun,t03d,opts);
%tval=fzero(@fun_strain,t03d,opts);
tval_strain=fzero(@fun_strain,t0_strain,opts);

if nargin>1
    figure(fig);
    plot(log10(dt),log10(left),'r--');hold on;
    plot(log10(dt),log10(right),'b--');
    plot(log10(dt),log10(left-right),'g--');
    plot(log10(tval),log10(rhs(tval)),'ks');
    axis([-10 -4 -8 8]);
    hold off;
    figure(fig+1);
    plot(log10(dt),log10(left_strain),'r--');hold on;
    plot(log10(dt),log10(right_strain),'b--');
    plot(log10(dt),log10(left_strain-right_strain),'g--');
    plot(log10(tval_1d),log10(rhs_strain(tval_strain)),'ks');
    axis([-10 -4 -8 8]);
    hold off;
end

function [tval tval_strain]=wave_solution_theta(t0,fig)

dt=logspace(-4,0,300); %row vector

%3d dipole excitation
left=ones(1,length(dt)); %row vector
right=rhs_theta(dt);
if nargin<1
    t0=[0.01 0.01];
end
if length(t0)<2
    t03d=t0;t01d=t0;
    t0_strain=t0;
else
    t03d=t0(1)*.85;t01d=t0(2)*.85;
    t0_strain=t0(2)*.85;
end
left_strain=ones(1,length(dt));
right_strain=rhs_strain_theta(dt);
opts = optimset('Display','off');
%tval=fzero(@fun_theta,[1e-6 0.1],opts);
tval=fzero(@fun_theta,t03d,opts);
%tval=fzero(@fun_strain_theta,t03d,opts);
%tval_strain=fzero(@fun_strain_theta,[1e-6 0.1],opts);
tval_strain=fzero(@fun_strain_theta,t0_strain,opts);

if nargin>1
    figure(fig);
    plot(log10(dt),log10(left),'r--');hold on;
    plot(log10(dt),log10(right),'b--');
    plot(log10(dt),log10(left-right),'g--');
    plot(log10(tval),log10(rhs_theta(tval)),'ks');
    axis(log10(2.71)*[-10 -4 -8 8]);
    hold off;
    figure(fig+1);
    plot(log10(dt),log10(left_strain),'r--');hold on;
    plot(log10(dt),log10(right_strain),'b--');
    plot(log10(dt),log10(left_strain-right_strain),'g--');
    plot(log10(tval_strain),log10(rhs_strain_theta(tval_strain)),'ks');
    axis(log10(2.71)*[-10 -4 -8 8]);
    hold off;
end

function result=fun(x)
result=rhs(x)-1;

function result=fun_theta(x)
result=rhs_theta(x)-1;

function result=fun_strain(x)
result=rhs_strain(x)-1;

function result=fun_strain_theta(x)
result=rhs_strain_theta(x)-1;

function result=rhs(x)
global E;global N;global Gamma;global nu;global r;
global modtype;
D = E*(1-nu)/(Gamma*(1+nu)*(1-2*nu));
y = 1./(4*D*sqrt(x.^2)); %hopefully will fix complaints
coeff = E.*y.^(7/2)/(r*(1-2*nu)*Gamma*pi^(3/2));
res=0;
for nx=1:N
    res=res+coeff.*exp(-nx.*y).*nx^(-3/2); %for ny=0 case
    for ny=1:N %symmetric under exchange ny <-> -ny gives factor of 2
        res=res+2*coeff.*exp(-(nx^2+ny^2).*y/nx).*(nx^2-ny^2)./nx^(7/2);
    end
end
fac=transmission_factor(E,modtype,nu);
result=res.*fac;

function result=rhs_theta(x)
global E;global N;global Gamma;global nu;global r;
global modtype;
if sum(isnan(x))>0
    result=zeros(size(x));
    result(isnan(x))=NaN;return;
end
D = E*(1-nu)/(Gamma*(1+nu)*(1-2*nu));
y = 1./(4*D*x); %hopefully will fix complaints
coeff = (1+nu)/(4*pi*(1-nu)*r); %simplified from E/(Gamma(1-2nu)*D);
res=0;tau=1e-1; %contraction time
for nx=1:N
    trQ=0;projQ=2;a2=nx.*y;a=sqrt(a2);
    res=res+coeff.*(W1(a,a2,trQ)+W2(a,a2,projQ));
    %need to add component that turns them off
    offidx=nx.*x - tau > 0; %which myocytes are turned off
    a2_d=nx.^2./(4*D*(nx.*x - tau));a_d=sqrt(a2_d);
    res=res-coeff.*(W1(a_d,a2_d,trQ,offidx)+W2(a_d,a2_d,projQ,offidx));
    for ny=1:N %symmetric under exchange ny <-> -ny gives factor of 2
        a2 = (nx+(ny^2/nx)).*y;a=sqrt(a2);
        trQ=0;projQ=(2*nx.^2-ny.^2)/(nx.^2+ny.^2);
        res=res+coeff.*(W1(a,a2,trQ)+W2(a,a2,projQ));
        %add off component via manual heaviside
        a2_d=(nx.^2+ny.^2)./(4*D*(nx.*x - tau));a_d=sqrt(a2_d);
        res=res-coeff.*(W1(a_d,a2_d,trQ,offidx)+W2(a_d,a2_d,projQ,offidx));
    end
end
fac=transmission_factor(E,modtype,nu);
result=res.*fac;

function result=rhs_strain(x)
global E;global N;global Gamma;global nu;global r;
global modtype;
D = E*(1-nu)/(Gamma*(1+nu)*(1-2*nu));
y = 1./(4*D*x); %hopefully will fix complaints
coeff = y.^(7/2)/(r*Gamma*pi^(3/2));
res=0;
for nx=1:N
    res=res+coeff.*exp(-nx.*y).*nx^(-3/2); %for ny=0 case
    for ny=1:N %symmetric under exchange ny <-> -ny gives factor of 2
        res=res+2*coeff.*exp(-(nx^2+ny^2).*y/nx).*(nx^2-ny^2)./nx^(7/2);
    end
end
fac=transmission_factor(E,modtype,nu);
result=res.*fac;

function result=rhs_strain_theta(x)
global E;global N;global Gamma;global nu;global r;
global modtype;
if sum(isnan(x))>0
    result=zeros(size(x));
    result(isnan(x))=NaN;return;
end
D = E*(1-nu)/(Gamma*(1+nu)*(1-2*nu));
y = 1./(4*D*x); %hopefully will fix complaints
coeff = 1./(4*pi*Gamma.*D*r); %simplified from E/(Gamma(1-2nu)*D);
res=0;tau=1e-1; %contraction time
for nx=1:N
    trQ=0;projQ=2;a2=nx.*y;a=sqrt(a2);
    res=res+coeff.*(W1(a,a2,trQ)+W2(a,a2,projQ));
    %need to add component that turns them off
    offidx=nx.*x - tau > 0; %which myocytes are turned off
    a2_d=nx.^2./(4*D*(nx.*x - tau));a_d=sqrt(a2_d);
    res=res-coeff.*(W1(a_d,a2_d,trQ,offidx)+W2(a_d,a2_d,projQ,offidx));
    for ny=1:N %symmetric under exchange ny <-> -ny gives factor of 2
        a2 = (nx+(ny^2/nx)).*y;a=sqrt(a2);
        trQ=0;projQ=(2*nx.^2-ny.^2)/(nx.^2+ny.^2);
        res=res+coeff.*(W1(a,a2,trQ)+W2(a,a2,projQ));
        %add off component via manual heaviside
        a2_d=(nx.^2+ny.^2)./(4*D*(nx.*x - tau));a_d=sqrt(a2_d);
        res=res-coeff.*(W1(a_d,a2_d,trQ,offidx)+W2(a_d,a2_d,projQ,offidx));
    end
end
fac=transmission_factor(E,modtype,nu);
result=res.*fac;

function result=fun1d(x)
result=rhs_dipole_1d(x)-lhs_dipole_1d(x);

function res=rhs_dipole_1d(x)
global E;global N;global Gamma;
En=E/Gamma;
y=x;
if N==0
    %polylog solution
else
    res=0;
    for n=1:N
        res=res+y.^(-3/2).*exp(-n./(4.*En.*y)).*(1./(2.*En.*y.*sqrt(n))-1./(n.^(3/2)));
    end
end

function res=lhs_dipole_1d(~)
global E;global r;global modtype;global Gamma;
fac=transmission_factor(E,modtype);
En=E/Gamma;
%res=r.*sqrt(16*pi*En)./(F(En,modtype).^2);
res=r.*sqrt(16*pi*En)./fac;

function fac=transmission_factor(E,type,nu)
%have to do something with the strain modulation
global Gamma;
if nargin<3
    nu=0;
end
if type==0
    fac=1;
elseif type==1
    fac=(0.75*L(E)).^2;
elseif type==2
    %fac=0.75.^2*E/((1-2*nu)*Gamma)*F(E).^2;
    %fac=0.75^2*F(E).^2;
    fac=0.75*F(E);
    %fac=0.75*E/((1-2*nu)*Gamma)*F(E);
elseif type==3
    %fac=E/(1-2*nu)/Gamma*M(E).^2;
    fac=M(E).^2;
    %fac=M(E);
end

function mod=M(E,n)
global E0;
if nargin<2
    n=6;
end
e=E/E0;
mod= 0.2 + .55./(1 + e.^n);

function mod=L(E,n)
global E0;

if nargin<2
    %n=1;
    %n=1.782606;
    %n=2.5;
    %n=3.5;
    n=4;
    %n=5;
end

e=E/E0;
mod=(4*e./((1+e).^2)).^n; %lorenzian fit

function mod=F(E)
global E0;
x=E/E0;
mod = tanh_log_hill(x);

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

function [val err]=expt_cellgel_vals(bool,norm)
if nargin<2
    norm=1/0.05;
end
if bool==0
    val=0;err=0;
else
    val=[0.3 0.84 2.6 10 40;norm*[0.022342 .049 .03649 .016 .0077]];
    dev=norm*[.017 .01009 .009 .0078 .012];
    err=[dev;-dev];
end

function result=tanh_log_hill(x,p)
if nargin<2
    p=1;
end
result=p*(0.5*tanh(-2*log(x)+1)+0.5).*(0.5*tanh(log(x)+3)+0.5);
%result=p*(0.5*tanh(-2*log(x)+1)+0.5); %flat profile
