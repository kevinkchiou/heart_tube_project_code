function [ output_args ] = plotConductionInterfere( fig )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    fig=1;
end

path(path,'./exptscripts'); %for some data scripts

adultNoSerumBPM=adultNoSerumHeart();
%adultSerumBPM=adultSerumHeart();
E4NoSerum=embryonicE4NoSerumHeart();
E6NoSerumBPM=embryonicE6NoSerumHeart();

%estimate bpm through average peak distances in autocorrelation fcn.
E4NoSerumBPM=zeros(3,2);temp=zeros(6,1);
for i=1:6
    t=E4NoSerum(:,1);
    datavec=E4NoSerum(:,1+i);
    temp(i)=estimateBPM(t,datavec);
end
%E4meanBPM=mean(temp);
for i=1:3
    E4NoSerumBPM(i,1)=mean([temp(i),temp(i+3)]);
    E4NoSerumBPM(i,2)=std([temp(i),temp(i+3)],0);
end
ETime(1:3,1)=[0,20,60]';

figure(fig);
%plot(adultNoSerumBPM(2:end,1)-20,adultNoSerumBPM(2:end,2),'k.-');hold on;%adult control
errorbar(adultNoSerumBPM(:,1)-20,adultNoSerumBPM(:,4),adultNoSerumBPM(:,5),'ro-');hold on;
%plot(adultSerumBPM(:,1)-20,adultSerumBPM(:,2),'rs-');
errorbar(ETime,E4NoSerumBPM(:,1),E4NoSerumBPM(:,2),'bo-');
errorbar(ETime,E6NoSerumBPM(:,1),E6NoSerumBPM(:,2),'go-'); %same time points between E4 and E6
hold off;
%axis([-20 60 0 300]);
axis([-25 65 0 220]);
xlabel('Time after +BGA (min)');ylabel('beats / min (BPM)');
%legend('Adult 25\mum BGA, -serum','Adult 100\mum BGA, +serum','E4 100\mum BGA, -serum','E6 100\mum BGA, -serum','Location','NorthEast');
%legend('Adult Control','Adult 25\mum BGA','E4 100\mum BGA','E6 100\mum BGA','Location','SouthEast');
legend('Adult 25\mum BGA','E4 100\mum BGA','E6 100\mum BGA','Location','NorthEast');
legend boxoff;
plot2svg('./conduct_interfere_BPM.svg',fig);
end

function out=estimateBPM(tvec,arvec)
acvec=autocorr(arvec,100);

sz=size(acvec);
%ensure a row column
if(sz(1)==1)
    %do nothing
elseif(sz(2)==1)
    acvec=acvec';
else
    error('fuck me');
end
%look beyond the first point
vec=acvec(2:end-1);
upvec=acvec(3:end);downvec=acvec(1:end-2);

%test element by element to find local maxima but boundaries don't count
maxima_idx=find(vec>upvec & vec>downvec);
minima_idx=find(vec<downvec & vec<upvec);

%find mean number of time points between maxima and minima
dn=mean([diff(maxima_idx),diff(minima_idx)]);
dt=dn*mean(diff(tvec)); %this is the average period in seconds.
out=60/dt; %convert to bpm
end
