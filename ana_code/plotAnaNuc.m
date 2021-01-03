%% Figures for thesis chapter on analytical solution
if 1
    baseSize=24;
    load('/Users/hhunt1/Documents/Nucleus/CutNuc1DData/cell7-1.mat')
    cellit=1;
    fits=dataStruct(cellit).fstfit{1};
    crtpts=dataStruct(cellit).crtpts{1};
    dataplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/data_plots/';
    anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
end
%% colours
cy = [240,228,66]/255;
co = [213 94 0]/255;
cb = [0 114 178]/255;
cg = [0 158 115]/255;
cp = [204 121 167]/255;
cgold= [230 159 0]/255;
csky = [86,180,233]/255;
%% Base computations
if 1
D=0.0174;
% D=0.001;
strD=num2str(D);
strD(strD=='.')=[];
timeInt=10;%(30*dataStruct(cellit).timeStep);
timeEnd=(1e4*dataStruct(cellit).timeStep);
time=0:timeInt:timeEnd;
[numSol,~,~]=nonlinearDiff(D,[],fits,crtpts,timeInt,timeEnd);
[anaSol1D,rfun,tfun]=sim_nucleus_cyl_spline(D,2,crtpts,fits,timeInt,timeEnd,0);
% [anaSol3Dline,h,anaSol3D,rfun,zfun,tfun]=sim_nucleus_cyl_wz_spline(D,20,2,10,crtpts,fits,timeInt,timeEnd,0);
[anaSol3D,~,~,~]=sim_nucleus_cyl_wz_spline_no_cp_c(D,20,2,{fits{1,1}},timeInt,timeEnd);
end
%% Compare numerical and analytical simulations
% Load cell7-1.mat

t1=[anaSol1D{:}];
t2=[anaSol3D{:}];
if 0
    figure
    surf(numSol','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    figure
    surf(t1,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
end
fpos=[400 250 660 300];
figure('pos',fpos)
hold on
p1=plot(time,t1(1,:),'LineWidth',4,'Color',cb);
p2=plot(time,t1(end,:),'LineWidth',2,'Color',co);
p3=plot(time,t2(1,501:2001:end),':','LineWidth',3,'Color',cy);
p4=plot(time,numSol(:,1),'--','LineWidth',2,'Color',cgold);
set(gca,'FontSize',0.7*baseSize)
ylabel('[Ca^{2+}] (\muM)','FontSize',baseSize*0.8)
xlabel('t (ms)','FontSize',baseSize*0.8)
legend([p2,p1,p3,p4],{'q_t(t)','1D analytical solution in the centre of the nucleus',...
    '3D analytical solution in the centre of the nucleus',...
    'numerical solution in the centre of the nucleus'},'FontSize',baseSize*0.6,...
    'Box','off')
saveas(gcf,strcat(anaplotsF,'compare_sol_D-',strD,'.eps'),'epsc')
if 1
% fpos=[400 250 560 340];
figure('pos',fpos)
plot(time,t1(1,:)-numSol(:,1)','LineWidth',2)
set(gca,'FontSize',0.7*baseSize)
title('Difference between 1D analytical and numerical','FontSize',baseSize*0.7,'FontWeight','normal')
ylabel('\Delta (\muM)','FontSize',baseSize*0.7)
xlabel('time (ms)','FontSize',baseSize*0.8)
axis([0 time(end) -10e-4 5e-4])
saveas(gcf,strcat(anaplotsF,'1D_diff_D-',strD,'.eps'),'epsc')

% fpos=[400 250 560 340];
figure('pos',fpos)
plot(time,t1(1,:)-t2(1,501:2001:end),'LineWidth',2)
set(gca,'FontSize',0.7*baseSize)
title('Difference between 1D and 3D analytical solutions (\muM)','FontSize',baseSize*0.7,'Fontweight','normal')
ylabel('\Delta (\muM)','FontSize',baseSize*0.7)
xlabel('time (ms)','FontSize',baseSize*0.8)
axis([0 time(end) -10e-4 5e-4])
saveas(gcf,strcat(anaplotsF,'ana_diff_D-',strD,'.eps'),'epsc')
end

%% Plot cross-section of nuc for 1D and 3D models at 180ms(?)
% load('cell7-1')
% cellit=1;
% pos=1;
% fits=dataStruct(cellit).fstfit{pos};
% crtpts=dataStruct(cellit).crtpts{pos};
anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
baseSize=38;
% tint=20;
% tend=5000;
a=2;l=20;
r=0:0.01:a;
z=0:0.01:l;
strD=strD;
strD(strD=='.')=[];
midNuc = cellfun(@(x)x(1,round(end/2)),anaSol3D,'Un',0);
midNuc = [midNuc{:}];
edgeNuc = cellfun(@(x)x(1,1),anaSol3D,'Un',0);
edgeNuc = [edgeNuc{:}];
[~,maxpt]=max(midNuc);
% tpts=round((mod(maxpt,200):200:timeEnd)/timeInt);
% tpts=27:5:45;
tpts=[27,32,42,52,62];
for tpt=tpts
    fpos=[400 250 560 340];
    figure('pos',fpos)
    surf([-r(end:-1:2),r],z,anaSol3D{tpt}([end:-1:1 2:(end)],:)','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    set(gca,'FontSize',baseSize*0.8)
    c=colorbar;
    c.Label.String='[Ca^{2+}] (\muM)';
    c.Label.FontSize=baseSize*0.8;
    cpos=c.Label.Position;
%     c.Label.Position=[cpos(1)-0.2 cpos(2:end)];
    title(strcat('time =',{' '},num2str(round((tpt-27)*timeInt)),{' '},'ms'),...
        'FontSize',baseSize,'FontWeight','normal')
    ylabel('z (\mu m)','FontSize',baseSize)
    xlabel('r (\mu m)','FontSize',baseSize)
    caxis([0.1 1])
    axis('tight')
    ax=gca;
    ax.Position=[0.14 0.27 0.6 0.6];
    saveas(gcf,strcat(anaplotsF,'nuc_crosssec_3D_D-',strD,'_',num2str(round((tpt-27)*timeInt)),...
        '_a-',num2str(a),'_l-',num2str(l),'p.eps'),'epsc')

    fpos=[400 250 560 340];
    figure('pos',fpos)
    surf([-r(end:-1:2),r],z,anaSol1D{tpt}([end:-1:1 2:(end)],ones(size(z)))','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    set(gca,'FontSize',baseSize*0.8)
    c=colorbar;
    c.Label.String='[Ca^{2+}] (\muM)';
    c.Label.FontSize=baseSize*0.8;
    cpos=c.Label.Position;
%     c.Label.Position=[cpos(1)-0.2 cpos(2:end)];
    title(strcat('time =',{' '},num2str(round((tpt-27)*timeInt)),{' '},'ms'),...
        'FontSize',baseSize,'FontWeight','normal')
    ylabel('z (\mu m)','FontSize',baseSize)
    xlabel('r (\mu m)','FontSize',baseSize)
    caxis([0.1 1])
    axis('tight')
    ax=gca;
    ax.Position=[0.14 0.27 0.6 0.6];
    saveas(gcf,strcat(anaplotsF,'nuc_crosssec_1D_D-',strD,'_',num2str(round((tpt-27)*timeInt)),...
        '_a-',num2str(a),'_l-',num2str(l),'p.eps'),'epsc')
end
%% Compare 1D and 3D simulations within areas of decreasing distance from the edge of the nucleus
anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
baseSize=24;
% tint=20;
% tend=5000;
time=0:timeInt:timeEnd;
a=2;l=20;
r=0:0.01:a;
z=0:0.01:l;
t1=[anaSol1D{:}];
edgeDist=1:50:201;

colours=[cg;cgold;cb;co;cp];
figure
lvals=cell(1,4);
for ii=edgeDist
box1D=cellfun(@(x)mean(x([ii:-1:1 2:ii])),anaSol1D,'Un',0);
box1D=[box1D{:}];
box3D=cellfun(@(x)mean(mean(x([ii:-1:1 2:ii],round(end/2+((-1000.5+ii*200.1/201):(1000.5-ii*200.1/201)))))),anaSol3D,'Un',0);
box3D=[box3D{:}];

hold on
plot(time,(box1D-box3D),'LineWidth',2);%,'Color',colours(ceil(ii/50),:));
% p2=plot(time,t1(end,:),'LineWidth',2,'Color',co);
% p3=plot(time,box3D,':','LineWidth',1,'Color',colours(ceil(ii/50),:));
% legend([p2,p1,p3],{'q(t)','1D solution','3D solution'},'FontSize',baseSize*0.6)
lvals{ceil(ii/50)}=char(strcat('Av. taken',{' '},num2str(round((ii-1)/100,2)),{' '},'\mum from edge'));
% if ii==1
%     [~,ptime]=max(box3D);
% end
end
% plot([time(ptime)*ones(1,2)],[-0.05,0.01],'--','LineWidth',2)
set(gca,'FontSize',baseSize*0.7)
ylabel('[Ca^{2+}] (\muM)','FontSize',baseSize*0.8)
xlabel('time (ms)','FontSize',baseSize*0.8)
title(strcat('Difference in av. [Ca^{2+}] of the 3D vs 1D model'),'FontSize',baseSize*0.9,'FontWeight','normal')
axis('tight')
legend(lvals,'FontSize',baseSize*0.7,'Box','off','Location','southeast')
saveas(gcf,strcat(anaplotsF,'1D_3D_edge_diff.eps'),'epsc')
%% Compare percentage difference between 1D and 3D simulations within areas of decreasing distance from the edge of the nucleus
anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
baseSize=24;
% tint=20;
% tend=5000;
time=0:timeInt:timeEnd;
a=2;l=20;
r=0:0.01:a;
z=0:0.01:l;
t1=[anaSol1D{:}];
edgeDist=1:50:201;

colours=[cg;cgold;cb;co;cp];
fpos=[400 250 760 300];
figure('pos',fpos)
lvals=cell(1,4);
hold on
for ii=edgeDist
box1D=cellfun(@(x)mean(x([ii:-1:1 2:ii])),anaSol1D,'Un',0);
box1D=[box1D{:}];
box3D=cellfun(@(x)mean(mean(x([ii:-1:1 2:ii],round(end/2+((-1000.5+ii*200.1/201):(1000.5-ii*200.1/201)))))),anaSol3D,'Un',0);
box3D=[box3D{:}];
plot(time,100*(box1D-box3D)./box1D,'LineWidth',2);%,'Color',colours(ceil(ii/50),:));
% p2=plot(time,t1(end,:),'LineWidth',2,'Color',co);
% p3=plot(time,box3D,':','LineWidth',1,'Color',colours(ceil(ii/50),:));
% legend([p2,p1,p3],{'q(t)','1D solution','3D solution'},'FontSize',baseSize*0.6)
lvals{ceil(ii/50)}=char(strcat('Av. taken',{' '},num2str(round((ii-1)/100,2)),{' '},'\mum from edge'));
end
set(gca,'FontSize',baseSize*0.7)
ylabel('\Delta[Ca^{2+}] (%)','FontSize',baseSize*0.8)
xlabel('t (ms)','FontSize',baseSize*0.8)
title(strcat('Difference in av. [Ca^{2+}] of the 1D vs 3D model'),'FontSize',baseSize*0.8,'FontWeight','normal')
legend(lvals,'FontSize',baseSize*0.7,'Location','southeast','AutoUpdate','off','Box','off')
plot([320,320],[-14,2],'--','LineWidth',2)
axis('tight')
saveas(gcf,strcat(anaplotsF,'1D_3D_edge_diff_pc.eps'),'epsc')
%% Plot caconc at 290ms 
%% Find difference in observed line-scan at various angles
% sliceAngles=[0,5,8,12,90];
sliceAngles=0:90;
l=10; a=2;
cellSlices=cell(size(sliceAngles));
hs=cellSlices;
slices=cellSlices;
sliceTTP=zeros(size(sliceAngles));
sliceFDHM=sliceTTP;
sliceMax=sliceTTP;
for cutangle=1:size(sliceAngles,2)
    [cellSlices{cutangle},hs{cutangle}]=cutnuc(sliceAngles(cutangle)*pi/180,l,a,anaSol3D,time);
    tempslice=cellSlices{cutangle};
    slices{cutangle}=[tempslice{:}];
    slicefit=fit(time',mean(slices{cutangle},1)','smoothingspline');
    newtime=time(1):0.01:time(end);
    mslice=slicefit(newtime);
    [slicemax,tempsliceTTP]=max(mslice);
    slicemin=min(mslice);
    sliceTTP(cutangle)=newtime(tempsliceTTP);
    sliceMax(cutangle)=slicemax;
    sliceFDHM(cutangle)=range(find(mslice>(slicemax+slicemin)/2))*0.01;
%     figure
%     surf(slices{cutangle},'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
%     view(90,90)
%     title(num2str(sliceAngles(cutangle)))

end

%% Plot difference in av. nuc transient at various angles
figure
hold on
numlines=5;
legstr=cell(1,numlines);
% if size(sliceAngles,2)>numlines+1
%     chooseSlice=1:floor(size(sliceAngles,2)/numlines):size(sliceAngles,2);
% end
chooseSlice=[0,5,8,12,90]+1;
for cutangle=1:numlines%chooseSlice
    plot(time,mean(slices{chooseSlice(cutangle)},1),'LineWidth',2)
%     legstr{ceil(cutangle/floor(size(sliceAngles,2)/numlines))}=strcat(num2str(sliceAngles(cutangle)),char(176));
    legstr{cutangle}=strcat(num2str(sliceAngles(chooseSlice(cutangle))),char(176));
end
set(gca,'FontSize',16)
leg=legend(legstr,'box','off');
title(leg,'Angle of line-scan','FontSize',baseSize*0.8)
title('Av. nuclear calcium transient','FontSize',baseSize,'FontWeight','normal')
xlabel('time (ms)','FontSize',baseSize*0.8)
ylabel('[Ca^{2+}] (\muM)','FontSize',baseSize*0.8)
axis('tight')
saveas(gcf,strcat(dataplotsF,'Av_trans_diff_D-',strD,'.eps'),'epsc')
%% Plot difference in med. nuc transient at various angles
figure
hold on
numlines=5;
legstr=cell(1,numlines+1);
if size(sliceAngles,2)>numlines+1
    chooseSlice=1:floor(size(sliceAngles,2)/numlines):size(sliceAngles,2);
end
for cutangle=chooseSlice
    plot(time,median(slices{cutangle},1),'LineWidth',2)
    legstr{ceil(cutangle/floor(size(sliceAngles,2)/numlines))}=strcat(num2str(sliceAngles(cutangle)),char(176));
end
set(gca,'FontSize',16)
leg=legend(legstr,'box','off');
title(leg,'Angle of line-scan','FontSize',baseSize*0.8)
title('Median nuclear calcium transient','FontSize',baseSize)
xlabel('time (ms)','FontSize',baseSize*0.8)
ylabel('[Ca^{2+}] (\muM)','FontSize',baseSize*0.8)
axis('tight')
%% Plot TTP vs angle
figure
plot(sliceAngles,sliceTTP,'LineWidth',2)
set(gca,'FontSize',16)
title(strcat('Av. nuclear calcium transient with D=',strD,'\mum^2/ms'),'FontSize',0.9*baseSize)
xlabel('angle (°)','FontSize',baseSize*0.8)
ylabel('TTP (ms)','FontSize',baseSize*0.8)
axis('tight')
saveas(gcf,strcat(dataplotsF,'TTPvAngle_D-',strD,'.eps'),'epsc')

%% Plot %TTP vs angle
figure
plot(sliceAngles,100*sliceTTP/sliceTTP(1),'LineWidth',2)
set(gca,'FontSize',16)
title(strcat('Av. nuclear calcium transient with D=',strD,'\mum^2/ms'),'FontSize',0.9*baseSize)
xlabel('angle (°)','FontSize',baseSize*0.8)
ylabel('TTP (% of TTP at 0°)','FontSize',baseSize*0.8)
% axis([0 90 0 100])
saveas(gcf,strcat(dataplotsF,'TTPvAngle_D-',strD,'pc.eps'),'epsc')
%% Plot %TTP vs angle scaled
figure
plot(sliceAngles,100*sliceTTP/sliceTTP(1),'LineWidth',2)
set(gca,'FontSize',baseSize*0.7)
title(strcat('Av. nuclear calcium transient with D=',num2str(D),'\mum^2/ms'),'FontSize',0.9*baseSize,'FontWeight','normal')
xlabel('angle (°)','FontSize',baseSize*0.8)
ylabel('TTP (% of TTP at 0°)','FontSize',baseSize*0.8)
axis([0 90 0 100])
saveas(gcf,strcat(dataplotsF,'TTPvAngle_D-',strD,'pcs.eps'),'epsc')
%% Plot FDHM vs angle
figure
plot(sliceAngles,sliceFDHM,'LineWidth',2)
set(gca,'FontSize',16)
title(strcat('Av. nuclear calcium transient with D=',strD,'\mum^2/ms'),'FontSize',0.9*baseSize)
xlabel('angle (°)','FontSize',baseSize*0.8)
ylabel('FDHM (ms)','FontSize',baseSize*0.8)
set(gca, 'YTickMode','manual')
set(gca, 'YTickLabel',num2str(get(gca,'YTick')'))
saveas(gcf,strcat(dataplotsF,'FDHMvAngle_D-',strD,'.eps'),'epsc')
%% Plot %FDHM vs angle
figure
plot(sliceAngles,100*sliceFDHM/sliceFDHM(1),'LineWidth',2)
set(gca,'FontSize',16)
title(strcat('Av. nuclear calcium transient with D=',strD,'\mum^2/ms'),'FontSize',0.9*baseSize)
xlabel('angle (°)','FontSize',baseSize*0.8)
ylabel('% of FDHM at 0°','FontSize',baseSize*0.8)
% set(gca, 'YTickMode','manual')
% set(gca, 'YTickLabel',num2str(get(gca,'YTick')'))
% axis([0 90 0 102])
saveas(gcf,strcat(dataplotsF,'FDHMvAngle_D-',strD,'pc.eps'),'epsc')
%% Plot %FDHM vs angle scaled
figure
plot(sliceAngles,100*sliceFDHM/sliceFDHM(1),'LineWidth',2)
set(gca,'FontSize',baseSize*0.7)
title(strcat('Av. nuclear calcium transient with D=',num2str(D),'\mum^2/ms'),'FontSize',0.9*baseSize,'FontWeight','normal')
xlabel('angle (°)','FontSize',baseSize*0.8)
ylabel('% of FDHM at 0°','FontSize',baseSize*0.8)
% set(gca, 'YTickMode','manual')
% set(gca, 'YTickLabel',num2str(get(gca,'YTick')'))
axis([0 90 0 110])
saveas(gcf,strcat(dataplotsF,'FDHMvAngle_D-',strD,'pcs.eps'),'epsc') 

%% 3: Effect of NE vs Dc & Effect of a vs Dc
if 1
anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
%load('/Users/hhunt1/Documents/huia/ana3DPS_tc1-6-11_midfixed')
saveBool=1;
baseSize=38;
numVars=5;
result_labels={'D_c','a','l','inclength','incheight','Max [Ca^{2+}]_{nuc} (\muM)','TTP (ms)','FDHM (ms)','Time to base (ms)','Max [Ca^{2+}]_{nuc av.} (\muM)','Av. TTP (ms)','Av. FDHM (ms)','Av. time to base (ms)'};
result_units={'\mum^2/ms','\mum','\mum','%','%','\muM','ms','ms','ms','\muM','ms','ms','ms'};
save_labels={'D_c','a','l','inclength','incheight','Max','ttpeak','FDHM','ttbase','MaxAv','ttpeakAv','FDHMAv','ttbaseAv'};
numpervar=round(exp(log(size(results,1))/numVars));
un=cell(1,numVars);
for unvar=1:numVars
    un{unvar}=unique(results(:,unvar));
end
resrangepreloop=find(results(:,1)==un{1}(6)); % Pref 1
% resrangepreloop=intersect(find(results(:,3)==un{3}(2)),resrangepreloop);
fixvars=[4,5];
fixvals=[1,7,5,6,6];
for fixvar=fixvars
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(fixvals(fixvar))),resrangepreloop);
end
for qual=[6:9]%,10:12]
    %[2,3]
    ii=2;
    jj=3;
    figure('position',[0,200,600,500])
    if 1
    try
        dc=reshape(results(resrangepreloop,qual),[numpervar,numpervar]);
    catch
        warning('wrong resrange')
    end
    x=unique(results(resrangepreloop,ii));
    y=unique(results(resrangepreloop,jj));
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    set(gca,'FontSize',baseSize*0.8)
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
    ax=gca;
    ax.Position=[0.22 0.25 0.6 0.6];
%     set(gca,'XScale','log')
    else
    x = results(resrangepreloop,ii);
    y = results(resrangepreloop,jj);
    dc= results(resrangepreloop,qual);
%     x = log(x);
    scatter(x,y,450,dc,'filled')
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    axis([min(x) max(x) min(y) max(y)])
%     set(gca,'XScale','log')
    end
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',baseSize)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize)
    c=colorbar;
    if qual==6 | qual== 10
    c.Ruler.TickLabelFormat = '%.2f';
    cpos=c.Label.Position;
    c.Label.Position=[cpos(1)-1 cpos(2:end)];
    c.Label.FontSize=baseSize*0.7;
    else
    c.Label.FontSize=baseSize*0.8;
    end
    c.Label.String=strcat(result_labels{qual});
    if saveBool
    saveas(gcf,strcat(anaplotsF,'hm-',save_labels{qual},'_',num2str(ii),...
        '_',num2str(jj)','_',...
        num2str(fixvars(1)),'-',num2str((results(resrangepreloop(1),fixvars(1)))),'_',...
        num2str(fixvars(2)),'-',num2str((results(resrangepreloop(1),fixvars(2)))),'c.eps'),'epsc')
    end
end
end

%% Plot cross-section of nuc over time
% load('cell7-1')
% fits=dataStruct(cellit).fstfit{pos};
% crtpts=dataStruct(cellit).crtpts{pos};
if 1
a=2;l=20;
r=0:0.01:a;
z=0:0.01:l;
tint=20;
tend=5000;
[sol,~,~,~]=sim_nucleus_cyl_wz_spline_no_cp_c(0.0174,l,a,{fits{1,1}},tint,tend);
end
midNuc = cellfun(@(x)x(1,round(end/2)),sol,'Un',0);
midNuc = [midNuc{:}];
edgeNuc = cellfun(@(x)x(1,1),sol,'Un',0);
edgeNuc = [edgeNuc{:}];
% [~,maxpt]=max(midNuc);
% tpts=round((mod(maxpt,200):200:tend)/tint);
baseSize=38;
% close all
for tpt=11:2:21%tpts(1:5)
    fpos=[400 250 560 340];
    figure('pos',fpos)
    surf([-r(end:-1:2),r],z,sol{tpt}([end:-1:1 2:(end)],:)','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    set(gca,'FontSize',baseSize*0.8)
    c=colorbar;
    c.Label.String = '[Ca^{2+}] (\muM)';
    c.Label.FontSize = baseSize*0.8;
    c.Ruler.TickLabelFormat = '%.2f';
    cpos=c.Label.Position;
    c.Label.Position=[cpos(1)-0.4 cpos(2:end)];
    title(strcat('time =',{' '},num2str((tpt-10)*tint),{' '},'ms'),'FontSize',baseSize,'FontWeight','normal')
    ylabel('z (\mu m)','FontSize',baseSize)
    xlabel('r (\mu m)','FontSize',baseSize)
    axis('tight')
    ax=gca;
    ax.Position=[0.14 0.27 0.6 0.6];
    saveas(gcf,strcat(anaplotsF,'nuc_crosssec_',num2str((tpt-10)*tint),...
        '_a-',num2str(a),'_l-',num2str(l),'.eps'),'epsc')
end
for tpt=11:2:21%tpts(1:5)
    fpos=[400 250 560 340];
    figure('pos',fpos)
    surf([-r(end:-1:2),r],z,sol{tpt}([end:-1:1 2:(end)],:)','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    set(gca,'FontSize',baseSize*0.8)
    c=colorbar;
    c.Label.String = '[Ca^{2+}] (\muM)';
    c.Label.FontSize = baseSize*0.8;
    c.Ruler.TickLabelFormat = '%.1f';
    title(strcat('time =',{' '},num2str((tpt-10)*tint),{' '},'ms'),'FontSize',baseSize,'FontWeight','normal')
    ylabel('z (\mu m)','FontSize',baseSize)
    xlabel('r (\mu m)','FontSize',baseSize)
%     caxis([0.7 1])
    caxis([0.1 0.9])
    axis('tight')
    ax=gca;
    ax.Position=[0.14 0.27 0.6 0.6];
    saveas(gcf,strcat(anaplotsF,'nuc_crosssec_',num2str((tpt-10)*tint),...
        '_a-',num2str(a),'_l-',num2str(l),'sc.eps'),'epsc')
end

%% Plot cross-section of nuc for different sizes at 180ms after transient
% load('cell7-1')
% cellit=1;
% pos=1;
% fits=dataStruct(cellit).fstfit{pos};
% crtpts=dataStruct(cellit).crtpts{pos};
anaplotsF='/Users/hhunt1/Documents/Nucleus/CutNuc1DData/ana_plots/';
baseSize=24;
avals=[2,2.4,2.4,1.6,1.6];
lvals=[20,24,16,16,24];
% hvals=[340,383.0667,383.0667,296.9333,296.9333];
hvals=[167.7333,176.3467,176.3467,159.12,159.12];
wvals=[560,658.7309,470.8131,470.8131,658.7309];

tpt=19;
tint=20;
tend=5000;
for ii=1:5
    a=avals(ii);
    l=lvals(ii);
r=0:0.01:a;
z=0:0.01:l;
[sol,~,~,~]=sim_nucleus_cyl_wz_spline_no_cp_c(0.0174,l,a,{fits{1,1}},tint,tend);
midNuc = cellfun(@(x)x(1,round(end/2)),sol,'Un',0);
midNuc = [midNuc{:}];
edgeNuc = cellfun(@(x)x(1,1),sol,'Un',0);
edgeNuc = [edgeNuc{:}];

    fpos=[400 250 wvals(ii) hvals(ii)];
    figure('pos',fpos)
    surf([-r(end:-1:2),r],z,sol{tpt}([end:-1:1 2:(end)],:)','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    hcb=colorbar;
    set(gca,'FontSize',baseSize*0.7)
    title(strcat(num2str(a),'x',num2str(l),'\mum'),'FontSize',baseSize*0.8,'FontWeight','normal')
    ylabel('z (\mu m)','FontSize',baseSize*0.8)
    xlabel('r (\mu m)','FontSize',baseSize*0.8)
    caxis([0.7 1])
    set(get(hcb,'Title'),'String','[Ca^{2+}] (\muM)','FontSize',baseSize*0.8)
    axis('tight')
    saveas(gcf,strcat(anaplotsF,'nuc_crosssec_2',num2str((tpt-10)*tint),...
        '_a-',num2str(a),'_l-',num2str(l),'p.eps'),'epsc')
end
