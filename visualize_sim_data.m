function [output]=visualize_sim_data(FILELIST)

output=0;
if nargin<1
    include=sprintf('grep ''output'' ');
    includediag=sprintf('grep ''diag'' ');
    exclude=sprintf('grep -v ''FILELIST\\|dat'' ');
    cmd=sprintf('ls | %s | %s > FILELIST.dat',include,exclude);
    cmddiag=sprintf('ls | %s | %s > FILELIST_DIAG.dat',includediag,exclude);
    system(cmd);system(cmddiag);
    FILELIST='FILELIST.dat';
    FILELISTDIAG='FILELIST_DIAG.dat';
end
%file with list of appropriate files
filelist=importdata(FILELIST);
flistdiag=importdata(FILELISTDIAG);
for i=1:length(filelist)
    FILENAME=char(filelist(i));
    if(~isempty(flistdiag))
        FILENAMEDIAG=char(flistdiag(i));
        myostate=plotdat(FILENAME,FILENAMEDIAG);
    else
        myostate=plotdat(FILENAME);
    end
    if i>1
        xprev=xnext;
        xnext=sum(myostate==0);
        v(i-1)=xprev-xnext;
    else
        xnext=sum(myostate==0);
    end
    M(i)=getframe;
    hold off;
end
vf=v(v>0);
v_avg=mean(vf)
movie(M,1);

end

function [myostate gridloc gridval]=plotdat(FILENAME,FILENAMEDIAG)

[gridval gridloc myostate]=readdat(FILENAME);
plot(1:length(gridval),gridval,'r.');hold on;
if nargin>1
    strain=readdiag(FILENAMEDIAG);
    plot(1:length(strain),strain,'b.');
end
refract=find(myostate==2);fire=find(myostate==1);%prime=find(myostate==0);
plot(gridloc(refract),zeros(length(gridloc(refract)),1),'ks');
plot(gridloc(fire),zeros(length(gridloc(fire)),1),'rs');
%plot(gridloc(prime),zeros(length(gridloc(prime)),1),'bs');
%axis equal;
end

function [gridval gridloc myostate]=readdat(FILENAME)

fp=fopen(FILENAME);
nums=fscanf(fp,'%d %d');

gridval=zeros(nums(1),1);
gridloc=zeros(nums(2),1);
myostate=zeros(nums(2),1);
for i=1:nums(1)
    gridval(i)=fscanf(fp,'%f',1);
end
for i=1:nums(2)
    gridloc(i)=fscanf(fp,'%d',1);
end
for i=1:nums(2)
    myostate(i)=fscanf(fp,'%d',1);
end
fclose(fp);

end

function gridval=readdiag(FILENAME)

fp=fopen(FILENAME);
nums=fscanf(fp,'%d');
gridval=zeros(nums(1),1);
for i=1:nums(1)
    gridval(i)=fscanf(fp,'%f',1);
end
fclose(fp);
end