%% H Hunt 2020
%% Extract results of D estimations from data
% load('/Users/hhunt1/Documents/Nucleus/cutNuc1DData/all_data_red.mat')
Dttp=[dataStruct(:).Dvalsm];
Dfdhm=[dataStruct(:).Dvalsfdhm];
nttp=[dataStruct(:).nttpeak];
nfdhm=[dataStruct(:).nfdhm];
cttp=[dataStruct(:).ttpeak];
cfdhm=[dataStruct(:).fdhm];
baseSize=24;
%% Stack estimates from ttpeak and FDHM
[~,Dbins]=histcounts([Dttp,Dfdhm],25);
DbinAv=(Dbins(1:(end-1))+Dbins(2:end))/2;
ttpcount=histcounts(Dttp,Dbins);
fdhmcount=histcounts(Dfdhm,Dbins);
figure
bar(DbinAv,[ttpcount;fdhmcount]','stacked')
legend({'Est. from TTP','Est. from FDHM'},'Box','off')
set(gca,'FontSize',16)
% title('Estimates of D','FontSize',baseSize)
xlabel('Estimate of D (\mum^2/ms)','FontSize',baseSize*0.8)
ylabel('Count','FontSize',baseSize*0.8)
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/Destimates_bar.eps','epsc')
%% Plot data on nuc dimensions
if 0
    load('/Users/hhunt1/Documents/Nucleus/cutNuc1DData/nucleus_dimensions.mat');
    baseSize=24;
end
% Plot hist of long axis
figure
histogram(nucDimensions(:,1),10)
set(gca,'FontSize',16)
title('Major axis of nucleus','FontSize',baseSize)
xlabel('Length (\mum)','FontSize',baseSize*0.8)
ylabel('Count','FontSize',baseSize*0.8)
axis('tight')
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/majoraxis.eps','epsc')
% Plot hist of two smaller axes
figure
hold on
histogram(nucDimensions(:,2),10)
histogram(nucDimensions(:,3),10)
set(gca,'FontSize',16)
title('Minor axes of nucleus','FontSize',baseSize)
xlabel('Length (\mum)','FontSize',baseSize*0.8)
ylabel('Count','FontSize',baseSize*0.8)
axis('tight')
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/minoraxes.eps','epsc')
%% Plot histogram of nuc and cyt FDHM
figure
hold on
hrange=(100*round(min([cfdhm,nfdhm]/100))):50:(100*round(max([cfdhm,nfdhm]/100)));
h1=histogram(nfdhm,hrange);
h2=histogram(cfdhm,hrange);
set(gca,'FontSize',16)
% title('FDHM','FontSize',baseSize,'FontWeight','normal')
xlabel('FDHM (ms)','FontSize',baseSize*0.8)
ylabel('frequency','FontSize',baseSize*0.8)
legend([h2,h1],{'cytosol','nucleus'},'FontSize',baseSize*0.7,'Box','off')
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/fdhm_hist.eps','epsc')
%% Plot histogram of nuc and cyt TTP
figure
hold on
hrange=(100*round(min([cttp,nttp]/100))):30:(100*round(max([cttp,nttp]/100)));
h1=histogram(nttp,hrange);
h2=histogram(cttp,hrange);
set(gca,'FontSize',16)
% title('TTP','FontSize',baseSize,'FontWeight','normal')
xlabel('TTP (ms)','FontSize',baseSize*0.8)
ylabel('frequency','FontSize',baseSize*0.8)
legend([h2,h1],{'cytosol','nucleus'},'FontSize',baseSize*0.7,'Box','off')
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/ttp_hist.eps','epsc')
%% Plot scatter of nuc vs cyt FDHM
figure
hold on
plot(cfdhm,nfdhm,'x','MarkerSize',10,'LineWidth',2)
plot([0 1200],[00, 1200],'LineWidth',2)
set(gca,'FontSize',16)
% title('FDHM','FontSize',baseSize,'FontWeight','normal')
xlabel('cytosol: FDHM (ms)','FontSize',baseSize*0.8)
ylabel('nucleus: FDHM (ms)','FontSize',baseSize*0.8)
axis([00 1200 0 1200])
legend({'data','y=x'},'FontSize',baseSize*0.6,'Location','southeast','Box','off')
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/fdhm_scatter.eps','epsc')

%% Plot scatter of nuc vs cyt TTP
figure
hold on
plot(cttp,nttp,'x','MarkerSize',10,'LineWidth',2)
plot([0 600],[00, 600],'LineWidth',2)
set(gca,'FontSize',16)
% title('FDHM','FontSize',baseSize,'FontWeight','normal')
xlabel('cytosol: TTP (ms)','FontSize',baseSize*0.8)
ylabel('nucleus: TTP (ms)','FontSize',baseSize*0.8)
axis([00 600 0 600])
legend({'data','y=x'},'FontSize',baseSize*0.6,'Location','southeast','Box','off')
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/ttp_scatter.eps','epsc')
%% Plot scatter of ttp vs fdhm D
figure
hold on
plot(Dttp,Dfdhm,'x','MarkerSize',10,'LineWidth',2)
plot([0.0 0.09],[0.0, 0.09],'LineWidth',2)
set(gca,'FontSize',16)
xlabel('D estimated from TTP (\mum^2/ms)','FontSize',baseSize*0.8)
ylabel('D estimated from FDHM (\mum^2/ms)','FontSize',baseSize*0.8)
legend({'Estimated D values','y=x'},'FontSize',baseSize*0.6,'Location','northwest','Box','off')
axis([0 0.07 0 0.07])
saveas(gcf,'/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/Dcalc_scatter.eps','epsc')

%% Plot each nucleus heatmap
sdir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/';
nextra=10;
for ii=1:size(dataStruct,1)
    tstep=dataStruct(ii).timeStep;
    nstep=dataStruct(ii).spaceStep;
    img=dataStruct(ii).smoothedImg;
    xrng=dataStruct(ii).xrng;
    trng=dataStruct(ii).init+(1:4000);
    nuc=dataStruct(ii).nucleus;
    nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
    cangle=abs(mod(dataStruct(ii).cutangle*180/pi,360));
    if cangle>180
        cangle=360-cangle;
    end
    figure('position',[0,200,600,500])
    surf((trng-trng(1))*tstep,(nrng-nrng(1))*nstep,img(nrng,trng),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(0,90)
    axis('tight')
    set(gca,'FontSize',0.7*baseSize)
    ylabel('position (\mum)','FontSize',baseSize*0.8)
    xlabel('time (ms)','FontSize',baseSize*0.8)
    ax=gca;
    ax.Position=[0.22 0.25 0.6 0.6];
    set(gca,'Xlim',[0 1000])
    title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
    c=colorbar;
    c.Ruler.TickLabelFormat = '%.1f';
    c.Label.String = 'F/F_0';
    c.Label.FontSize = baseSize*0.8;
    cpos=c.Label.Position;
    c.Label.Position=[cpos(1)-0.6 cpos(2:end)];
%     axlims=get(gcf,'CurrentAxes');
%     newx=round(axlims.YLim(1)):round(axlims.YLim(end));
%     set(gca,'YTick',newx(1:round(end/3):end));
%     xpos=get(xl,'Position');
%     set(xl,'Position',[xpos(1)+2 xpos(2)+2e2 xpos(3:end)])
    saveas(gcf,strcat(sdir,'cell',num2str(ii),'_nuc1_extra',num2str(nextra),'.eps'),'epsc')
    if ~isempty(dataStruct(ii).nucleus2)
        nuc=dataStruct(ii).nucleus2;
        nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
        figure('position',[0,200,600,500])
        surf((trng-trng(1))*tstep,(nrng-nrng(1))*nstep,img(nrng,trng),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
        view(0,90)
        axis('tight')
        set(gca,'FontSize',0.7*baseSize)
        ylabel('position (\mum)','FontSize',baseSize*0.8)
        xlabel('time (ms)','FontSize',baseSize*0.8)
        ax=gca;
        ax.Position=[0.22 0.25 0.6 0.6];
        set(gca,'Xlim',[0 1000])
        title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
        c=colorbar;
        c.Ruler.TickLabelFormat = '%.1f';
        c.Label.String = 'F/F_0';
        c.Label.FontSize = baseSize*0.8;
        cpos=c.Label.Position;
        c.Label.Position=[cpos(1)-0.6 cpos(2:end)];
        saveas(gcf,strcat(sdir,'cell',num2str(ii),'_nuc2_extra',num2str(nextra),'.eps'),'epsc')
    end
end
%% Plot each nucleus heatmap long
sdir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/';
nextra=0;
baseSize=48;
for ii=1%1:size(dataStruct,1)
    tstep=dataStruct(ii).timeStep;
    nstep=dataStruct(ii).spaceStep;
    img=dataStruct(ii).smoothedImg;
    xrng=dataStruct(ii).xrng;
    trng=dataStruct(ii).init+(1:4000);
    nuc=dataStruct(ii).nucleus;
    nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
    cangle=abs(mod(dataStruct(ii).cutangle*180/pi,360));
    if cangle>180
        cangle=360-cangle;
    end
    figure('position',[400 250 760 340])
    surf((trng-trng(1))*tstep,(nrng-nrng(1))*nstep,img(nrng,trng),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    axis('tight')
    set(gca,'Xlim',[0 1000])
    view(0,90)
    c=colorbar;
    set(gca,'FontSize',0.7*baseSize)
    ylabel('position (\mum)','FontSize',baseSize*0.8)
    xlabel('time (ms)','FontSize',baseSize*0.8)
%     ax=gca;
%     ax.Position=[0.22 0.25 0.6 0.6];
%     title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
    cpos=c.Position;
    c.Position=[cpos(1) cpos(2:end)];
    c.FontSize=0.8*baseSize;
    set(get(c,'Title'),'String','      (\muM)','FontSize',baseSize*0.7)
    clpos=c.Title.Position;
    c.Title.Position=[clpos(1)+20 clpos(2)-30 clpos(3:end)];
%     c.Label.String = 'F/F_0';
%     c.Label.FontSize = baseSize*0.8;
%     cpos=c.Label.Position;
%     c.Label.Position=[cpos(1)-1.7 cpos(2:end)];
%     axlims=get(gcf,'CurrentAxes');
%     newx=round(axlims.YLim(1)):round(axlims.YLim(end));
%     set(gca,'YTick',newx(1:round(end/3):end));
%     xpos=get(xl,'Position');
%     set(xl,'Position',[xpos(1)+2 xpos(2)+2e2 xpos(3:end)])
    saveas(gcf,strcat(sdir,'cell_long',num2str(ii),'_nuc1_extra',num2str(nextra),'.eps'),'epsc')
end
%% Plot each nucleus mesh
sdir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/';
nextra=0;
for ii=1%1:size(dataStruct,1)
    tstep=dataStruct(ii).timeStep;
    nstep=dataStruct(ii).spaceStep;
    img=dataStruct(ii).smoothedImg;
    xrng=dataStruct(ii).xrng;
    trng=dataStruct(ii).init+(1:3000);
    nuc=dataStruct(ii).nucleus;
    nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
    cangle=abs(mod(dataStruct(ii).cutangle*180/pi,360));
    if cangle>180
        cangle=360-cangle;
    end
    figure
    surf((trng(1:50:end)-trng(1))*tstep,(nrng(1:end)-nrng(1))*nstep,img(nrng(1:end),trng(1:50:end)),'FaceColor','interp','FaceLighting','gouraud');
    view(-45,50)
    set(gca,'FontSize',0.7*baseSize)
    ylabel('\mum','FontSize',baseSize*0.8)
    xlabel('time (ms)','FontSize',baseSize*0.8)
%     title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
    c=colorbar;
    c.Label.String = 'F/F_0';
    c.Label.FontSize = baseSize*0.8;
    axlims=get(gcf,'CurrentAxes');
    newx=round(axlims.YLim(1)):round(axlims.YLim(end));
    set(gca,'YTick',newx(1:round(end/3):end));
    saveas(gcf,strcat(sdir,'mesh_cell',num2str(ii),'_nuc1_extra',num2str(nextra),'.eps'),'epsc')
    if ~isempty(dataStruct(ii).nucleus2)
        nuc=dataStruct(ii).nucleus2;
        nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
        figure
        surf((trng(1:50:end)-trng(1))*tstep,(nrng-nrng(1))*nstep,img(nrng,trng(1:50:end)),'FaceColor','interp','FaceLighting','gouraud');
        view(-45,50)
        set(gca,'FontSize',0.7*baseSize)
        ylabel('\mum','FontSize',baseSize*0.8)
        xlabel('time (ms)','FontSize',baseSize*0.8)
        title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
        c=colorbar;
        c.Label.String = 'F/F_0';
        c.Label.FontSize = baseSize*0.8;
        saveas(gcf,strcat(sdir,'mesh_cell',num2str(ii),'_nuc2_extra',num2str(nextra),'.eps'),'epsc')
    end
end

%% Plot each nucleus contour
sdir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/data_plots/';
nextra=10;
numcs=12;
for ii=[4,1]%1:size(dataStruct,1)
    tstep=dataStruct(ii).timeStep;
    nstep=dataStruct(ii).spaceStep;
    img=dataStruct(ii).smoothedImg;
    xrng=dataStruct(ii).xrng;
    trng=dataStruct(ii).init+(1:3000);
    nuc=dataStruct(ii).nucleus;
    nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
    cangle=abs(mod(dataStruct(ii).cutangle*180/pi,360));
    if cangle>180
        cangle=360-cangle;
    end
    figure
    contourf((trng-trng(1))*tstep,(nrng-nrng(1))*nstep,img(nrng,trng),numcs);
%     view(90,-90)
    set(gca,'FontSize',0.7*baseSize)
    ylabel('\mum','FontSize',baseSize*0.8)
    xlabel('time (ms)','FontSize',baseSize*0.8)
    title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
    c=colorbar;
    c.Label.String = 'F/F_0';
    c.Label.FontSize = baseSize*0.8;
    saveas(gcf,strcat(sdir,'cont_cell',num2str(ii),'_nuc1_extra',num2str(nextra),'.eps'),'epsc')
    if 0%~isempty(dataStruct(ii).nucleus2)
        nuc=dataStruct(ii).nucleus2;
        nrng=(xrng(nuc(1))-nextra):xrng(nuc(end))+nextra;
        figure
        contourf((trng-trng(1))*tstep,(nrng-nrng(1))*nstep,img(nrng,trng),numcs);
%         view(90,-90)
        set(gca,'FontSize',0.7*baseSize)
        ylabel('\mum','FontSize',baseSize*0.8)
        xlabel('time (ms)','FontSize',baseSize*0.8)
        title(strcat(num2str(cangle),'°'),'FontWeight','normal','FontSize',baseSize*0.8)
        c=colorbar;
        c.Label.String = 'F/F_0';
        c.Label.FontSize = baseSize*0.8;
%         saveas(gcf,strcat(sdir,'cont_cell',num2str(ii),'_nuc2_extra',num2str(nextra),'.eps'),'epsc')
    end
end

%% Find curve
%img(nrng,trng)
curvec=zeros(size(nrng));
for ii=1:size(nrng,2)
    tempvar=find(img(nrng(ii),trng)>1.4,1,'first');
    if ~isempty(tempvar)
        curvec(ii)=tempvar;
    end
end
curves=find(curvec(1:round(end/2))==0,1,'last')+1;
% curvee=find(curvec(round(end/2):end)==0,1,'first');

curvec=curvec(curves:end);%(curvee+round(end/2)));
% 

r=zeros(size(curvec));
for ii=2:(size(curvec,2)-1)
    tempd=cross([ii,curvec(ii),0]-[ii-1,curvec(ii-1),0],[ii+1,curvec(ii+1),0]-[ii-1,curvec(ii-1),0]);
    r(ii)=norm([ii-1,curvec(ii-1)])*norm([ii,curvec(ii)])*norm([ii+1,curvec(ii+1)])/2/norm(tempd);
end