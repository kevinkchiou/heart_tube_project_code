function [out_func] = engler_data_fit(meth)
%ENGLER_DATA_FIT() takes Engler's data and creates a fitted function to it
%via method meth.
%   Detailed explanation goes here
if nargin<1
    meth=1;
end

%Engler's raw data obtained via datathief in units of kpa and strain %
data=[1.0 11.8;8.1 15.1;11.0 12.9;17.1 2.9;23.2 3.3;34.2 3.8];
E0=10;

%normalize
data(:,1)=data(:,1)/E0;
data(:,2)=data(:,2)/max(data(:,2));
x=data(:,1)';y=data(:,2)';

if meth==1
    cs=spline([0 x 50/E0],[0 y 0]); %clamped spline
    cs.coefs
    xx=linspace(0,50/E0,100);
    plot(x,y,'bs',xx,ppval(cs,xx),'r-');
elseif meth==2
    params=[6 1];
    p=lsqcurvefit(@e_hill_fcn,params,x,y);
    xx=linspace(0,50/E0,100);
    plot(x,y,'bs',xx,e_hill_fcn(p,xx),'r-');
elseif meth==3
    p=1;
    p=lsqcurvefit(@tanh_log_hill,p,x,y);
    xx=linspace(0,50/E0,100);
    set(plot(x,y,'b.'),'LineWidth',2,'MarkerSize',30);hold on;
    set(plot(xx,tanh_log_hill(p,xx),'r-'),'LineWidth',2);
    xlabel('E / E^*');ylabel('Strain transmission');
    legend('Experiment','Fit Curve','Location','NorthEast');hold off;
end
end

function result=tanh_log_hill(p,x)
result=p(1)*(0.5*tanh(-2*log(x)+1)+0.5).*(0.5*tanh(log(x)+3)+0.5);
%result=p(1)*(0.5*tanh(-2*log(x)+1)+0.5);
end

function result=e_hill_fcn(params,xdata)
p=params;x=xdata;
result=0.2+0.8./(1+exp(p(1).*(x-p(2))));
end