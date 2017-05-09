function [ out ] = singleCellDataFit(fig,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    noplot=1;
end
if nargin<2
    opt=1; %use weighted function
end

%here, x, x0 are vector of parameters: [Es, Ec, eps_trans]
x0=[1 1 1]; %initially Es=Ec=E^*, eps_trans=1 is overall strain transmission scale of 1 (full eshelby)
xmin=[0.5 0.5 0.1];xmax=[2.0 2.0 1.0]; %maximum and minimum values for parameters
if opt==1
    out=lsqnonlin(@fun_error,x0,xmin,xmax);
elseif opt==0
    out=lsqnonlin(@fun_noerror,x0,xmin,xmax);
else
    error('Error: invalid option entry')
end

if noplot==0
    figure(fig);close(fig);figure(fig);
    %we can plot the resulting curve here with a full set of E's instead of
    %the reduced set just from Stephanie's values
end

end

function out=fun_error(x,opt)
if nargin<2
    opt=0;
end
%x is the vector of parameters, out is the output vector of values
Es=x(1);Ec=x(2);
if numel(x)>2
    eps_trans=x(3);
end
%Use Stephanie's values to fit
[ec ece]=expt_cellgel_vals(1);
E=ec(1,:);V = ec(2,:); W = max(ece);% Stephanie scaled values and error weights
out=zeros(size(ec));
%can't vectorize in substrate stiffness computation since each one is a
%simulation - run a (short) loop!
for i=1:numel(E)
    eps=callComsoleFcn(E(i),Es,Ec,eps_trans);
    if opt==1
        %use error-weighted values
        if W(i)>0
            out(i)=(eps-V(i))/W(i);
        else
            %case where weight has a funny value...
            W(i)=mean([W(1:i-1),W(i+1:end)]);
        end
        out(i)=(eps - V(i))/W(i);
    elseif opt==0
        out(i)=eps-V(i);
    else
        error('Error: undefined option in least squares function evaluation');
    end
end
end

function strain_tr=callComsoleFcn(E,Es,Ec,eps_trans)
strain_tr=JasonsComsoleFcn(E,Es,Ec,eps_trans);
end

%Stephanie's values
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
    err(1,val(2,:)-dev<0)=val(2,val(2,:)-dev<0); %prevent error bars < 0
end
end