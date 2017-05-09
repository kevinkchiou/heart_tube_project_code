function [ output_args ] = halfspace_singlecell_plots(fig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    fig=1;
end
figure(fig);

dat=csvread('../jason/halfspace_comsole/E_t_vs_s_fit.csv',1);
semilogx(dat(:,1),dat(:,2),'r-');hold on; %plot Jason's data

%plot the other data
[ec ece]=expt_cellgel_vals(1,1); %normalization of 1
errorbar(ec(1,:),ec(2,:),ece(1,:),ece(2,:),'ks');hold off;
xlabel(' E / E^* ');ylabel('Tr \epsilon');
axis([0.1,100,0,0.06]);
end