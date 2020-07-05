if 0
    baseSize=24;
    load('/Users/hhunt1/Documents/Nucleus/CutNuc1DData/cell7-1.mat')
    cellit=1;
    fits=dataStruct(cellit).fstfit{1};
    crtpts=dataStruct(cellit).crtpts{1};
    anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
end
%% For period = 1000 ms
inclength=1;
incheight=1;
% period=180;%3000;
periodl=1000;
dfit=dilatepp(fits{1,1},inclength,incheight);
nxvalsl=dfit.breaks;
nxvalsl=nxvalsl;
nxvalsl(nxvalsl>(nxvalsl(1)+periodl))=[];
nyvalsl=ppval(dfit,nxvalsl);
tempxvals=[nxvalsl(end-2),nxvalsl(end-1),nxvalsl(1)+periodl];
tempyvals=[nyvalsl(end-2),nyvalsl(end-1),nyvalsl(1)];
tempf=fit(tempxvals',tempyvals','smoothingspline','SmoothingParam',0.07);
nxvalsl=[nxvalsl,(nxvalsl(end-1)+nxvalsl(end))/2];%,3*(nxvals(end-1)+nxvals(end))/4];
nyvalsl=[nyvalsl,0];%,0];
nyvalsl((end-2):end)=tempf(nxvalsl((end-2):end));
nxvalsl=nxvalsl(1:(end-1));nyvalsl=nyvalsl(1:(end-1));
%% For period = 180 ms
inclength=1;
incheight=1;
% period=180;%3000;
periods=180;
dfit=dilatepp(fits{1,1},inclength,incheight);
nxvalss=dfit.breaks;
nxvalss=nxvalss;
nxvalss(nxvalss>(nxvalss(1)+periods))=[];
nyvalss=ppval(dfit,nxvalss);
tempxvals=[nxvalss(end-2),nxvalss(end-1),nxvalss(1)+periods];
tempyvals=[nyvalss(end-2),nyvalss(end-1),nyvalss(1)];
tempf=fit(tempxvals',tempyvals','smoothingspline','SmoothingParam',0.07);
nxvalss=[nxvalss,(nxvalss(end-1)+nxvalss(end))/2];%,3*(nxvalss(end-1)+nxvalss(end))/4];
nyvalss=[nyvalss,0];%,0];
nyvalss((end-2):end)=tempf(nxvalss((end-2):end));
%%

repNum=5;
repxvalsl=repmat(nxvalsl,[1,repNum])+reshape(repmat(0:periodl:periodl*(repNum-1),[size(nxvalsl,2),1]),[1,repNum*size(nxvalsl,2)]);
repxvalss=repmat(nxvalss,[1,repNum*2])+reshape(repmat(0:periods:periods*(repNum*2-1),[size(nxvalss,2),1]),[1,repNum*2*size(nxvalss,2)]);
repxvals=[repxvalsl(1:end),repxvalss+repxvalsl(end)];
repyvalsl=repmat(nyvalsl,[1,repNum]);
repyvalss=repmat(nyvalss,[1,repNum*2]);
repyvals=[repyvalsl(1:end),repyvalss];
af=fit(repxvals',repyvals','smoothingspline','SmoothingParam',0.07);
fitmult=coeffvalues(af);

if 1
    figure
    plot(af,repxvals,repyvals)
end
%%
timeInt=10;
timeEnd=repxvals(end);
D=0.0174;
% D=0.001;
[anaSol3Dmult,~,~,~]=sim_nucleus_cyl_wz_spline_no_cp_c(D,20,2,{fitmult},timeInt,timeEnd);

t3Dmult=[anaSol3Dmult{:}];
time=0:timeInt:timeEnd;
%%
fpos=[400 250 860 340];
    figure('pos',fpos)
hold on
% plot(af,repxvals,repyvals)
plot(repxvals(1):repxvals(end),ppval(fitmult,repxvals(1):repxvals(end)),'LineWidth',2)
p3=plot(time,t3Dmult(1,501:2001:end),'LineWidth',2);
axis([0,repNum*periodl+2*repNum*periods,0,incheight*1.05])
axis([0,6726,0,incheight*1.05])
set(gca,'FontSize',baseSize*0.7)
% title(strcat('Multiple transients at T=',num2str(periodl),{' '},'(ms)'),'FontSize',baseSize*)
ylabel('[Ca^{2+}] (\mu M)','FontSize',0.8*baseSize)
xlabel('time (ms)','FontSize',0.8*baseSize)
legend({'cytosol','nucleus'},'Location','south','FontSize',0.7*baseSize)
saveas(gcf,strcat(anaplotsF,'multi_D-00',num2str(D*1e4,'%05.f'),'_T-',num2str(periodl),...
     '_lh-',num2str(1/inclength*100-100),'-',num2str(incheight*100-100),'.eps'),'epsc')
saveas(gcf,strcat(anaplotsF,'multi_D-00',num2str(D*1e4,'%05.f'),'_T-',num2str(periodl),...
    '_lh-',num2str(1/inclength*100-100),'-',num2str(incheight*100-100),'.fig'))
%%
if 0
figure
hold on
end
skint=2001*round(periodl/timeInt);
strans=zeros(repNum,round(periodl/timeInt)+1);
for ii=1:(repNum-1)
    strans(ii,:)=t3Dmult(1,(skint*(ii-1)+501):2001:(skint*ii+501));
    if 0
    plot(0:timeInt:periodl,strans(ii,:));
    end
end
skint=2001*round(periods/timeInt);
strans2=zeros(2*repNum,round(periods/timeInt)+1);
for ii=(1:2*(repNum-1))
    strans2(ii,:)=t3Dmult(1,(skint*ii+501)+((skint*(ii-1)+501):2001:(skint*ii+501)));
    if 0
    plot(0:timeInt:periodl,strans(ii,:));
end
end
%%
figure
plot(0:timeInt:periodl,strans(2:(end-2),:)'-strans(1:(end-3),:)','LineWidth',2)
set(gca,'FontSize',baseSize*0.7)
title('Difference between transients','FontSize',baseSize*0.8,'FontWeight','normal')
ylabel('[Ca^{2+}] (\mu M)','FontSize',0.8*baseSize)
xlabel('time (ms)','FontSize',0.8*baseSize)
legend(strcat('transient',{' '},num2str((2:(repNum-1))'),'-',num2str((1:(repNum-2))')),'FontSize',0.7*baseSize,'Location','southeast')
saveas(gcf,strcat(anaplotsF,'multi_diff_D-00',num2str(D*1e4,'%05.f'),'_T-',num2str(periodl),...
     '_lh-',num2str(1/inclength*100-100),'-',num2str(incheight*100-100),'.eps'),'epsc')
 
 %%
 figure
 hold
plot(0:timeInt:(timeInt*((round(periodl/timeInt)+1)*(repNum-3)-1)),reshape(strans(2:(end-2),:)',[1,(round(periodl/timeInt)+1)*(repNum-3)])-reshape(strans(1:(end-3),:)',[1,(round(periodl/timeInt)+1)*(repNum-3)]),'LineWidth',2)
plot((timeInt*((round(periodl/timeInt)+1)*(repNum-3)-1))+...
    (0:timeInt:(timeInt*((round(periods/timeInt)+1)*(repNum-3)-1)))...
    ,reshape(strans2(2:(end-2),:)',[1,(round(periods/timeInt)+1)*(2*repNum-3)])-...
    reshape(strans2(1:(end-3),:)',[1,(round(periods/timeInt)+1)*(2*repNum-3)]),'LineWidth',2)
title('Difference between subsequent nuclear transients','FontSize',0.9*baseSize)
ylabel('[Ca^{2+}] (\mu M)','FontSize',0.8*baseSize)
xlabel('time (ms)','FontSize',0.8*baseSize)
axis('tight')
% saveas(gcf,strcat(anaplotsF,'multi_difflong_D-00',num2str(D*1e4,'%05.f'),'_T-',num2str(periodl),...
%      '_lh-',num2str(1/inclength*100-100),'-',num2str(incheight*100-100),'.eps'),'epsc')
 
 %% Measure FDHM of each trans
 [av,ap]=findpeaks(-ppval(fitmult,repxvals(1):repxvals(end)));
 lastTransTimePts=round(ap(end-1)/timeInt):round(ap(end)/timeInt);
 lastTransTime=(lastTransTimePts-lastTransTimePts(1))*timeInt;
 lastTrans=zeros(6,size(lastTransTimePts,2));
 temptrans=[anasoln{lastTransTimePts}];
 lastTrans(1,:)=temptrans(1,501:2001:end);
 lastTrans(4,:)=temptrans(end,501:2001:end);
 temptrans=[anasoll{lastTransTimePts}];
 lastTrans(2,:)=temptrans(1,501:2001:end);
 lastTrans(5,:)=temptrans(end,501:2001:end);
 temptrans=[anasolh{lastTransTimePts}];
 lastTrans(3,:)=temptrans(1,501:2001:end);
 lastTrans(6,:)=temptrans(end,501:2001:end);
 figure
 plot(lastTrans')
 
 newTrans=zeros(6,range(lastTransTime)+1);
 newTime=(lastTransTime(1):lastTransTime(end));
 for ii=1:6
    newFit=fit(lastTransTime',lastTrans(ii,:)','smoothingspline','SmoothingParam',0.07);
    newTrans(ii,:)=newFit(newTime);
 end
     
 maxv=zeros(1,6);
 minv=maxv;
 fdhm=maxv;
for ii=1:6
    [maxv(ii),~]=max(newTrans(ii,:));
    [minv(ii),~]=min(newTrans(ii,:));
    fdhm(ii)=range(find(newTrans(ii,:)>(maxv(ii)+minv(ii))/2));
end
figure
plot(newTime,newTrans')