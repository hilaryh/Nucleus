%% Plot numDiffPS results

% Gridded params
savedir='/Volumes/UniWD/AnaNucleus/numPS/';
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
    savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
end
pos=1;
numVars=6;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'log \mum^2/ms','log \mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
%%
for ii=1:(numVars-1)
%     unii=unii(2:end)';
    for jj=(ii+1):numVars
        otherVars=setdiff(1:6,[ii,jj]);
        % Make sure all other params are fixed
        % 1: D_c 2: D_b 3: b_t 4: K 5: NPC 6: a
        un1=unique(results(:,otherVars(1)))';
        un1=un1(2:end);
        un2=unique(results(:,otherVars(2)))';
        un2=un2(3);%un2(2:end);
        un3=unique(results(:,otherVars(3)))';
        un3=un3(3);%un3(2:end);
        un4=unique(results(:,otherVars(4)))';
        un4=un4(3);%un4(2:end);
%         unjj=unjj(3:end)';
        for l2=un2
            for l1=un1
                for l3=un3
                    for l4=un4
                        resrange1=find(results(:,otherVars(1))==l1);
                        resrange2=find(results(resrange1,otherVars(2))==l2);
                        resrange3=find(results(resrange2,otherVars(3))==l3);
                        resrange4=find(results(resrange3,otherVars(4))==l4);
                        resrange=resrange1(resrange2(resrange3(resrange4)));
                        for qual = numVars+(1:4)
                            x = results(resrange,ii);
                            y = results(resrange,jj);
                            dc= results(resrange,qual);
                            if ii==1 || ii==2
                                x = log(x);
                            end
                            if jj==2
                                y = log(y);
                            end

                            figure('position',[0,200,560,450])
                            scatter3(x,y,dc,1000,dc,'filled')
                            set(gca,'FontSize',16)
                            view(90,-90)
                            axis('square')
                            colorbar
                            title(strcat(result_labels{qual}),'FontSize',25)
                            ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',20)
                            xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',20)
                            saveas(gcf,strcat(savedir,'g-scatter_pos-',num2str(pos),'_',num2str(ii),...
                                '_',num2str(jj),'-',num2str(l1),'-',num2str(l2),'-',num2str(l3),'-',num2str(l4),'_',num2str(qual),'.eps'),'epsc')
        %                     saveas(gcf,strcat(savedir,'scatter_pos-',num2str(pos),'_',num2str(ii),'-',num2str(kk),...
        %                         '_',num2str(jj),'-',num2str(ll),'_',num2str(qual),'.fig'),'epsc')
                        end
                    end
                end
        close all
            end
        end
    end
end

%% Gridded params in 3d scatter plots 
% Gridded params
savedir='/Volumes/UniWD/AnaNucleus/numPS/';
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
    savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
end
pos=1;
numVars=6;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'log(\mum^2/ms)','log(\mum^2/ms)','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
for ii=1:(numVars-1)
%     unii=unii(2:end)';
    for jj=(ii+1):numVars
        otherVars=setdiff(1:6,[ii,jj]);
        % Make sure all other params are fixed
        % 1: D_c 2: D_b 3: b_t 4: K 5: NPC 6: a
        for kk=otherVars(otherVars>jj)
            fixVars=setdiff(otherVars,[kk]);
            un1=unique(results(:,fixVars(1)))';
            un1=un1(2);
            un2=unique(results(:,fixVars(2)))';
            un2=un2(2);%un2(2:end);
            un3=unique(results(:,fixVars(3)))';
            un3=un3(2);%un3(2:end);
            resrange1=find(results(:,fixVars(1))==un1);
            resrange2=find(results(resrange1,fixVars(2))==un2);
            resrange3=find(results(resrange1(resrange2),fixVars(3))==un3);
            resrange=resrange1(resrange2(resrange3));
            for qual = numVars+(1:4)
                x = results(resrange,ii);
                y = results(resrange,jj);
                z = results(resrange,kk);
                dc= results(resrange,qual);
                if ii==1 || ii==2
                    x = log(x);
                end
                if jj==2
                    y = log(y);
                end
                if kk==2 || kk==1
                    z = log(z);
                end

                figure('position',[0,200,560,450])
                scatter3(x,y,z,1000,dc,'filled')
                set(gca,'FontSize',16)
                view(-17,12)
                axis('square')
                title(strcat(result_labels{qual}),'FontSize',25)
                ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',20)
                xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',20)
                zlabel(strcat(result_labels(kk),{' '},'(',result_units(kk),')'),'FontSize',20)
                colorbar
                saveas(gcf,strcat(savedir,'g-scatter3_pos-',num2str(pos),'_',save_labels{qual},'_',num2str(ii),...
                    '_',num2str(jj),'_',num2str(kk),'_fix-',num2str(un1),'-',num2str(un2),'-',num2str(un3),'.eps'),'epsc')
%                     saveas(gcf,strcat(savedir,'scatter_pos-',num2str(pos),'_',num2str(ii),'-',num2str(kk),...
%                         '_',num2str(jj),'-',num2str(ll),'_',num2str(qual),'.fig'),'epsc')
                        end
        close all
                    end
                end
end
            
%% Make similar grid but with no extras
Dvals = unique(results(:,1));

%% Plot each parameter with sorted quals
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'log(\mum^2/ms)','log(\mum^2/ms)','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numa=4;
numd=5;
resrange=find(results(:,5)>0);
resrange=intersect(find(results(:,6)==4.5),resrange);
resrangeb=intersect(find(results(:,3)>1),resrange);
resranged{1}=intersect(resrangeb,find(results(:,1)<results(:,2))); 
resranged{2}=intersect(resrangeb,find(results(:,1)==results(:,2)));
resranged{3}=intersect(resrangeb,find(results(:,1)>results(:,2)));
resrangepreloop=intersect(find(results(:,3)<1),resrange);
for ii=7:9
    [sortqualb, iiorderb] = sort(results(resrangeb,ii));
    for kk=1:3
        [sortquald{kk}, iiorderd{kk}] = sort(results(resranged{kk},ii));
        rsiid{kk}=results(resranged{kk}(iiorderd{kk}),1:6);
        rsiid{kk}(:,1:2)=log(rsiid{kk}(:,1:2));
    end
    [sortqualnb, iiordernb] = sort(results(resrangepreloop,ii));
    rsiib=results(resrangeb(iiorderb),1:6);
    rsiinb=results(resrangepreloop(iiordernb),1:6);
    rsiib(:,1:2)=log(rsiib(:,1:2));
    rsiinb(:,1:2)=log(rsiinb(:,1:2));
    figure
    scatLabel={'no buf','D_b>D_c','D_b=D_c','D_b<D_c'};
    for jj=1:5
        subplot(numd,numa,4*jj-3)
%         scatter(sortqualb,rsiib(:,jj),'.','MarkerEdgeAlpha',0.1)
        hold on
        scatter(sortqualnb,rsiinb(:,jj),200,'.','MarkerEdgeAlpha',0.1)
        ylabel(result_labels{jj})
        if jj==1
            title(result_labels{ii})
        elseif jj==5
            xlabel(scatLabel{1})
        end
        for kk=1:3
            subplot(numd,numa,4*jj-3+kk)
            scatter(sortquald{kk},rsiid{kk}(:,jj),200,'.','MarkerEdgeAlpha',0.1)
            if jj==5
                xlabel(scatLabel{kk+1})
            end
        end
        
        end
%     title(result_labels{ii})
end

%% Figures for chapter
if 0
%% 2: Data vs Model
% Load cell7-1.mat
cellit=6;
mccpos=dataStruct(cellit).mccpos{1};
baseSize=24;
meanpeakDist=mean(diff(mccpos));
figure
% Data
% subplot(1,2,1)
thisPeak1=dataStruct(cellit).thisPeak{1}{1};
cyto=median(dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(dataStruct(cellit).nucleus(1:10)-15),...
                (max(round(mccpos(1)-0.1*meanpeakDist)+dataStruct(cellit).init,1):...
            min(round(mccpos(1)+0.9*meanpeakDist)+dataStruct(cellit).init,size(dataStruct(cellit).smoothedImg,2)))),1);
mccpos=dataStruct(cellit).mccpos{1};
meanpeakDist=mean(diff(mccpos));
nuco=median(dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(dataStruct(cellit).nucleus((end/2-3):(end/2+3))),...
                (max(round(mccpos(1)-0.1*meanpeakDist)+dataStruct(cellit).init,1):...
            min(round(mccpos(1)+0.9*meanpeakDist)+dataStruct(cellit).init,size(dataStruct(cellit).smoothedImg,2)))),1);
npeakMin=mean(nuco(end-10:end));
npeakMax=max(nuco);
global peakMin; global peakMax;
peakMin=mean(cyto(end-10:end));
peakMax=max(cyto);
% nuc=(0.626*nuco+0.1*peakMax-0.726*peakMin)/(peakMax-peakMin);
t1=(1:size(thisPeak1,2))*dataStruct(cellit).timeStep;
hold on
if 0
plot(t1,thisPeak1)
plot(t1,nuc)
end
plot(t1,nuco,'LineWidth',2)
plot(t1,cyto,'LineWidth',2)

set(gca,'FontSize',baseSize*0.7)
title('Data','FontSize',baseSize*0.9,'FontWeight','normal')
 xlabel('time (ms)','FontSize',baseSize*0.8)
 ylabel('[Ca^{2+}] (F/F_0)','FontSize',baseSize*0.8)
%% Simulation
% subplot(1,2,2)
figure
hold on
fits=dataStruct(cellit).fstfit{1};
crtpts=dataStruct(cellit).crtpts{1};
if 1
% [cutca,~,~]=nonlinearDiff(dataStruct(cellit).Dvalsm(1),[],fits,crtpts,dataStruct(cellit).timeStep,max(t1)-dataStruct(cellit).timeStep);
[cutca,~,~]=nonlinearDiff(0.0176,[],fits,crtpts,dataStruct(cellit).timeStep,max(t1)-dataStruct(cellit).timeStep);
plot(t1,cutca(:,1)+npeakMin-peakMin,'LineWidth',2)
plot(t1,q_spline(t1,fits,crtpts),'LineWidth',2)
% axis([0 max(t1) 0 1])
else
[cutca2,~,~]=nonlinearDiff(dataStruct(cellit).Dvalsm(1),[],fits,crtpts,30*dataStruct(cellit).timeStep,3000*dataStruct(cellit).timeStep);
plot(t1,cutca2(:,1)+npeakMin-peakMin,'LineWidth',2);
plot(t1,(q_spline(t1,fits,crtpts)),'LineWidth',2)%*(peakMax-peakMin)-0.1*peakMax+0.726*peakMin)/0.626);
end
 set(gca,'FontSize',baseSize*0.8)
legend('nucleus','cytosol','FontSize',baseSize,'Box','off')
 title('Simulation','FontSize',baseSize*0.9,'FontWeight','normal')
 xlabel('time (ms)','FontSize',baseSize*0.8)
end

%% 3: Effect of NE vs Dc & Effect of a vs Dc
% Load PS1
% load('/Users/hhunt1/Documents/huia/parPS1')
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
baseSize=42;
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
resrange=find(results(:,5)>0);
resrange=intersect(find(results(:,6)==un{6}(5)),resrange); % Pref 2
resrangepreloop=intersect(find(results(:,4)==0),resrange);
for fixvar=[2,3]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(1)),resrangepreloop);
end
if 0
    
    set(gca,'YScale','log')
    set(gca,'FontSize',sqSize*0.8)
    axlims=get(gcf,'CurrentAxes');
    newy=10.^(floor(log10(axlims.YLim(1))):floor(log10(axlims.YLim(end))));
    set(gca,'YTick',newy(1:round(end/3):end));
    ytickformat('%.0e');
	axis('tight')
    xl=xlabel(strcat(result_labels(23-lvar),{' '},'(',result_units(23-lvar),')'),'FontSize',sqSize*0.8);
    yl=ylabel(strcat(result_labels(1),{' '},'(',result_units(1),')'),'FontSize',sqSize*0.8);
    ypos=get(yl,'Position');
    set(yl,'Position',[ypos(1)+1e-2 ypos(2:end)])
    xpos=get(xl,'Position');
    set(xl,'Position',[xpos(1) xpos(2)+1e-5 xpos(3:end)])
    c=colorbar;
    c.Label.String=strcat(result_labels(qual),{' '},'(',result_units(qual),')');
    c.Label.FontSize=sqSize*0.8;
end
for qual=7:9
    ii=1;
    jj=5;
    
    figure('position',[0,200,600,500])
    if 1
    dc=reshape(results(resrangepreloop,qual),[numpervar,numpervar]);
    x=unique(results(resrangepreloop,ii));
    y=unique(results(resrangepreloop,jj));
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XScale','log')
    ax=gca;
    ax.Position=[0.22 0.25 0.6 0.6];
    else
    x = results(resrangepreloop,ii);
    y = results(resrangepreloop,jj);
    dc= results(resrangepreloop,qual);
    scatter(x,y,4500,dc,'filled')
    axis('square')
    axis([min(x-2) max(x+2) min(y-0.25) max(y+0.25)])
    set(gca,'XScale','log')
    end
    set(gca,'FontSize',baseSize*0.7)
    c=colorbar;
    c.Label.String=strcat(result_labels{qual},{' '},'(',result_units{qual},')');
    c.Label.FontSize=baseSize*0.8;
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',baseSize*0.8)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize*0.8)
    axlims=get(gcf,'CurrentAxes');
    newx=10.^(floor(log10(axlims.XLim(1))):floor(log10(axlims.XLim(end))));
    try
        set(gca,'XTick',newx(1:round(end/3):end));
    catch
        warning('Odd bug, set manually')
    end
    saveas(gcf,strcat(savedir,'chapterfig_NE_pos-',num2str(pos),'_',save_labels{qual},'_',num2str(ii),...
        '_',num2str(jj),'.eps'),'epsc')
end

resrange=find(results(:,5)==un{5}(6)); % Pref 1
resrange=intersect(find(results(:,6)>0),resrange);
resrangepreloop=intersect(find(results(:,4)==0),resrange);
for fixvar=[2,3]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(1)),resrangepreloop);
end
for qual=7:9
    ii=1;
    jj=6;
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
    set(gca,'FontSize',baseSize*0.7)
    ax=gca;
    ax.Position=[0.22 0.25 0.6 0.6];
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XScale','log')
    else
    x = results(resrangepreloop,ii);
    y = results(resrangepreloop,jj);
    dc= results(resrangepreloop,qual);
    x = log(x);
    scatter(x,y,4500,dc,'filled')
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    axis([min(x) max(x) min(y) max(y)])
    set(gca,'XScale','log')
    end
    c=colorbar;
    c.Label.String=strcat(result_labels{qual},{' '},'(',result_units{qual},')');
    c.Label.FontSize=baseSize*0.8;
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',baseSize*0.8)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize*0.8)
    axlims=get(gcf,'CurrentAxes');
    newx=10.^(floor(log10(axlims.XLim(1))):floor(log10(axlims.XLim(end))));
    try
        set(gca,'XTick',newx(1:round(end/3):end));
    catch
        warning('Odd bug. Set manually')
    end
    saveas(gcf,strcat(savedir,'chapterfig_a_pos-',num2str(pos),'_',save_labels{qual},'_',num2str(ii),...
        '_',num2str(jj),'.eps'),'epsc')
end
%% 6: nuclear transient generated for various values of D_c
Dc=logspace(-4,0,7);
figure
hold on
legText=result_labels;
time=0:2.6:4000;
inclVals=[1:4,6];
for ii=inclVals
%     [cutca{ii},~,~]=nonlinearDiff_hl(Dc(ii),[],1,1,fits,crtpts,2.6,4000);
    plot(time,cutca{ii}(:,1),'LineWidth',2)
    legText{ii}=strcat('D_c=',num2str(round(Dc(ii),1,'significant')),'\mum^2/ms');
end
set(gca,'FontSize',baseSize*0.8)
title('Simulated nuclear calcium','FontSize',baseSize*1.25)
ylabel(strcat('[Ca^{2+}]',{' '},'(',result_units(7),')'),'FontSize',baseSize)
xlabel(strcat('time',{' '},'(','ms',')'),'FontSize',baseSize)
legend(legText{inclVals},'FontSize',0.8*baseSize,'box','off','Location','best')

%% 8: Effect of buffers
% Load PS1
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=result_labels;
tText={'D_b<D_c','D_b>D_c'};
dval=4;
resrange=find(results(:,5)==un{5}(5));
resrange=intersect(find(results(:,6)==un{6}(5)),resrange);
resrangepreloop=intersect(find(results(:,1)==un{1}(dval)),resrange);
for qual=7:9
    figure
    for dbit=[1,2]
        resrange=intersect(find(results(:,2)==un{1}(2*dbit+5-dval)),resrangepreloop);
        ii=3;
        figure
        hold on
        count=0;
        maxP=0;minP=1e4;
        for jj=1:2:size(un{4},1)
            count=count+1;
            resrangeK{count}=intersect(find(results(:,4)==un{4}(jj)),resrange);
            legText{count}=strcat('K=',num2str(round(un{4}(jj),2,'significant')),'\muM');
            x = results(resrangeK{count},ii);
            dc= results(resrangeK{count},qual);
            if max(dc)>maxP
                maxP=max(dc);
            end
            if min(dc)<minP
                minP=min(dc);
            end
            if count<8
                plot(x,dc,'-o','MarkerSize',5,'LineWidth',2)
            else
                plot(x,dc,'-x','MarkerSize',5,'LineWidth',2)
            end
        end
        set(gca,'XScale','log')
        set(gca,'FontSize',baseSize*0.8)
        if dbit==1
            title(strcat(result_labels{qual},':',{' '},tText{dbit}),'FontSize',baseSize*1.25)
        ylabel(strcat(result_labels(qual),{' '},'(',result_units(qual),')'),'FontSize',baseSize)
        else
            title(strcat(tText{dbit}),'FontSize',baseSize*1.25)
        end
        axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
            minP maxP])
        xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize)
        if dbit==2
            if qual==10 
            legend(legText,'Location','best','Box','off','FontSize',baseSize*0.6)
        else
            legend(legText,'Location','best','Box','off','FontSize',baseSize*0.6)
            end
        end
%         saveas(gcf,strcat(savedir,'chapterfig4_pos-',num2str(pos),'D_c',num2str(un{1}(dval)),'_',save_labels{qual},'_',num2str(ii),...
%             '_',num2str(dbit),'.eps'),'epsc')
    end
end

%% 8: Effect of buffers - heatmaps
% Load PS1_bk12
savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
pos=1;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
un=cell(1,6);
baseSize=42;
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=result_labels;
tText={'D_b<D_c','D_b>D_c'};
dval=1;
numpervar=300;
% resrange=1:size(results,1);
resrange=find(results(:,5)==un{5}(5));
resrange=intersect(find(results(:,6)==un{6}(5)),resrange);
resrangepreloop=intersect(find(results(:,1)==un{1}(dval)),resrange);
for qual=7:9
    for dbit=[1,2]
        resrange=intersect(find(results(:,2)==un{1}(2*dbit+5-dval)),resrangepreloop);
        ii=3;
        jj=4;
        fpos=[400 250 600 500];
        figure('pos',fpos)
        hold on
        count=0;
        maxP=0;minP=1e4;
            x = unique(results(resrange,ii));
            y = unique(results(resrange,jj));
            dc= reshape(results(resrange,qual),[numpervar,numpervar]);
            surf(x,y,dc,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'FontSize',baseSize*0.7)
        ax=gca;
        ax.Position=[0.22 0.25 0.6 0.6];
        axlims=get(gcf,'CurrentAxes');
        newy=10.^(floor(log10(axlims.XLim(1))):floor(log10(axlims.XLim(end))));
        set(gca,'XTick',newy(1:round(end/3):end));
        xtickformat('%.0e');
        xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize*0.8)
        yl=ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',baseSize*0.8);
        ypos=get(yl,'Position');
        set(yl,'Position',[ypos(1)+1e-5 ypos(2:end)])
        c=colorbar;
        if qual==7
        c.Ruler.TickLabelFormat = '%.1f';
        end
        c.Label.String=strcat(result_labels(qual),{' '},'(',result_units(qual),')');
        c.Label.FontSize=baseSize*0.8;
        axis('tight')
        if qual==7
           lpos=c.Label.Position; 
           c.Label.Position=[lpos(1)-0.5,lpos(2),lpos(3)];
        end
        if dbit==1
            title(strcat(tText{dbit}),'FontSize',baseSize,'FontWeight','normal')
        else
            title(strcat(tText{dbit}),'FontSize',baseSize,'FontWeight','normal')
        end
        saveas(gcf,strcat(savedir,'chapterfig4_hm_pos-',num2str(pos),'D_c',num2str(un{1}(dval)),'_',save_labels{qual},'_',num2str(ii),...
            '_',num2str(dbit),'.eps'),'epsc')
    end
end
%% Buffer heatmaps
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=result_labels;
tText={'D_b<D_c','D_b>D_c'};
dval=4;
resrange=find(results(:,5)==un{5}(5));
resrange=intersect(find(results(:,6)==un{6}(5)),resrange);
resrangepreloop=intersect(find(results(:,1)==un{1}(dval)),resrange);
for qual=7:9
    for dbit=[1,2]
        resrange=intersect(find(results(:,2)==un{1}(2*dbit+5-dval)),resrangepreloop);
        ii=3;
        figure
        hold on
        count=0;
        maxP=0;minP=1e4;
            x = unique(results(resrange,ii));
            y = unique(results(resrange,4));
            dc= reshape(results(resrange,qual),[numpervar,numpervar]);
            if max(dc)>maxP
                maxP=max(dc);
            end
            if min(dc)<minP
                minP=min(dc);
            end
            surf(x,y,dc,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
            view(90,90)
            axis('square')
            axis('tight')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'FontSize',baseSize*0.8)
        colorbar
        title(strcat(result_labels{qual},{' '},'(',result_units{qual},')',':',{' '},tText{dbit}),'FontSize',baseSize,'FontWeight','normal')
        ylabel(strcat(result_labels(4),{' '},'(',result_units(4),')'),'FontSize',baseSize*0.8)
        xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize*0.8)
        saveas(gcf,strcat(savedir,'chapterfig4_pos_hm-',num2str(pos),'D_c',num2str(un{1}(dval)),'_',save_labels{qual},'_',num2str(ii),...
            '_',num2str(dbit),'.eps'),'epsc')
    end
end
%% d_c v D_b v fdhm
% Load PS1
savedir='/Volumes/UniWD/AnaNucleus/numPS/';
savebool=1;
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
    savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
end
pos=1;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
resrange=find(results(:,5)==un{5}(5));
resrangepreloop=intersect(find(results(:,6)==un{6}(5)),resrange);
% fixval=[1,1,11,8,1,1];
fixcam=[1,1,1,6,1,1];
fixne =[1,1,5,4,1,1];
fixcalretc=[1,1,9,9,1,1];
fixval=fixcam;
for fixvar=[3,4]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(fixval(fixvar))),resrangepreloop);
end
for qual=7:9
    ii=1;
    jj=2;
    
    figure('position',[0,200,560,450])
    dc=reshape(results(resrangepreloop,qual),[numpervar,numpervar]);
    x=unique(results(resrangepreloop,ii));
    y=unique(results(resrangepreloop,jj));
    psurf=0;
    if psurf
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
    else
    [C,h]=contourf(x,y,dc);
%     clabel(C,h,'Color','w','FontSize',baseSize*0.8)
    end
    set(gca,'FontSize',baseSize*0.8)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    axis([1e-4 1e2 1e-4 1])
%     axis('square')
    axis('tight')
    clabel(C,h,'manual','Color','w','FontSize',baseSize*0.8)
    title(strcat(result_labels{qual}),'FontSize',1.25*baseSize)
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',baseSize)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize)
    if psurf
        colorbar
        saveas(gcf,strcat(savedir,'chapterfig_Dc_Db_s_pos-',num2str(pos),'_',save_labels{qual},'_',num2str(ii),...
        '_',num2str(jj),'.eps'),'epsc')
    else
        if savebool
        saveas(gcf,strcat(savedir,'chapterfig_Dc_Db_pos-',num2str(pos),'_',save_labels{qual},'_',num2str(ii),...
        '_',num2str(jj),'c1.eps'),'epsc')
%         figure('position',[0,200,560,450])
%         [C,h]=contourf(dc');
%         clabel(C,h,'Color','w','FontSize',baseSize*0.8)
%         set(gca,'FontSize',16)
%         saveas(gcf,strcat(savedir,'chapterfig_Dc_Db_pos-',num2str(pos),'_',save_labels{qual},'_',num2str(ii),...
%         '_',num2str(jj),'c2.eps'),'epsc')
        end
    end
    
end



%% d_c v D_b v fdhm for every K
% Load PS1
savedir='/Volumes/UniWD/AnaNucleus/numPS/';
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
    savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
end
pos=1;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end

psurf=1;
resrange=find(results(:,5)==un{5}(5));
resrangepreloop=intersect(find(results(:,6)==un{6}(5)),resrange);
fixval=[1,1,2,8,1,1];
for fixvar=[3]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(fixval(fixvar))),resrangepreloop);
end
numa=4;numd=3;
for qual=7:9
    ii=1;
    jj=2;
    
    figure('units','normalized','outerposition',[0 0 1 1])
for kk=1:numpervar
    resrangeend=intersect(find(results(:,4)==un{4}(kk)),resrangepreloop);
    subplot(numd,numa,kk)
    dc=reshape(results(resrangeend,qual),[numpervar,numpervar]);
    x=unique(results(resrangeend,ii));
    y=unique(results(resrangeend,jj));
    if psurf
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
    else
    [C,h]=contourf(x,y,dc');
    end
    set(gca,'FontSize',16)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    title(strcat('K=10^{',num2str(-6+kk),'}\muM'),'FontSize',25)
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',20)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',20)
    c=colorbar;
    if qual==7
        caxis([0 1])
    elseif qual==8
        caxis([32 300])
    elseif qual==9
        caxis([280 2000])
    end
end
    if psurf
        c.Label.String=(strcat(result_labels{qual},{' '},'(',result_units{qual},')'));
    end
end
if psurf
%     saveas(gcf,strcat(savedir,'chapterfig_Dc_Db_compare_K_pos-',num2str(pos),...
%         '_',save_labels{qual},'b_t-',num2str(round(results(resrangeend(1),3),2)),'.eps'),'epsc')
end

%% Compare numerical and analytical simulations
% Load cell7-1.mat
cellit=1;
numpervar=round(exp(log(size(results,1))/6));
fits=dataStruct(cellit).fstfit{1};
crtpts=dataStruct(cellit).crtpts{1};
[cutca,~,~]=nonlinearDiff(1e-3,[],fits,crtpts,30*dataStruct(cellit).timeStep,1e4*dataStruct(cellit).timeStep);
[sol,rfun,tfun]=sim_nucleus_cyl_spline(1e-3,2,crtpts,fits,30*dataStruct(cellit).timeStep,1e4*dataStruct(cellit).timeStep,0);
t1=[sol{:}];
figure
surf(cutca','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
view(90,90)
figure
surf(t1,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
view(90,90)

%% Chapterfig 10-13: Effect of buffers Db-Dc grid
% Load PS1
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=cell(size(un{4}));
tText={'D_b<D_c','D_b>D_c'};
dval=3;
resrange=find(results(:,5)==un{5}(5));
resrange=intersect(find(results(:,6)==un{6}(5)),resrange);
resrangepreloop=(resrange);
for qual=8%7:9
    figure('units','normalized','outerposition',[0 0 1 1])
    for dbit=1:((numpervar-5)*(numpervar-6))
        resrange=intersect(find(results(:,1)==un{1}(ceil(dbit/(numpervar-5)))),resrangepreloop);
        resrange=intersect(find(results(:,2)==un{2}(mod(dbit-1,(numpervar-5))+1)),resrange);
        count=0;
        ii=3;
        subplot(numpervar-6,numpervar-5,dbit)
        hold on
        for jj=1:size(un{4},1)
            count=count+1;
            resrangeK{count}=intersect(find(results(:,4)==un{4}(jj)),resrange);
            legText{count}=strcat('K=',num2str(round(un{4}(jj),2,'significant')),'\muM');
            x = results(resrangeK{jj},ii);
            dc= results(resrangeK{jj},qual);
            if jj<8
                h=plot(x,dc,'-o','MarkerSize',5,'LineWidth',2);
            else
                plot(x,dc,'-x','MarkerSize',5,'LineWidth',2)
            end
        end
        set(gca,'XScale','log')
        title(strcat('D_c=',num2str(round(un{1}(ceil(dbit/(numpervar-5))),1,'significant')),'D_b=',num2str(round(un{2}(mod(dbit-1,numpervar-5)+1),1,'significant'))))
        if qual==9
            axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
            0 3600])
        elseif qual==8
            axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
            0 520])
        elseif qual==7
            axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
            0 1])
        end
        ylabel(strcat(result_labels(qual),{' '},'(',result_units(qual),')'))%,'FontSize',20)
        xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'))%,'FontSize',20)
        
% 
        numberOfXTicks = 4;
        xData = get(h,'XData');
        set(gca,'Xtick',logspace(log(xData(2))/log(10),log(xData(end))/log(10),numberOfXTicks))


    end
    if qual==7
        legend(legText,'Location','best','Box','off','FontSize',12)
    else
        legend(legText,'Location','best','Box','off','FontSize',12)
    end
%     saveas(gcf,strcat(savedir,'Dc_Db_K_b_t','_',save_labels{qual},...
%         '.eps'),'epsc')
end
%% Chapterfig 10-13: Effect of buffers Db-Dc grid - heatmaps
% Load PS1
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
% legText=cell(size(un{4}));
% tText={'D_b<D_c','D_b>D_c'};
dval=3;
resrange=find(results(:,5)==un{5}(5));
resrange=intersect(find(results(:,6)==un{6}(5)),resrange);
resrangepreloop=(resrange);
for qual=8%:9
    figure('units','normalized','outerposition',[0 0 1 1])
    for dbit=1:((numpervar-5)*(numpervar-6))
        resrange=intersect(find(results(:,1)==un{1}(ceil(dbit/(numpervar-5)))),resrangepreloop);
        resrange=intersect(find(results(:,2)==un{2}(mod(dbit-1,(numpervar-5))+1)),resrange);
        count=0;
        ii=3;
        jj=4;
        subplot(numpervar-6,numpervar-5,dbit)
        hold on
        x = unique(results(resrange,ii));
        y = unique(results(resrange,jj));
        dc= reshape(results(resrange,qual),[numpervar,numpervar]);
        h=surf(x,y,dc,'FaceColor','interp','FaceLighting','gouraud','EdgeColor','none');
        view(0,90)
        set(gca,'XScale','log')
        set(gca,'YScale','log')
%         title(strcat('D_c=',num2str(round(un{1}(ceil(dbit/(numpervar-5))),1,...
%             'significant')),{' '},'D_b=',num2str(round(un{2}(mod(dbit-1,numpervar-5)+1),...
%             1,'significant'))),'FontWeight','normal')
        ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'))%,'FontSize',20)
        xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'))%,'FontSize',20)
        c=colorbar;
        c.Label.String=strcat(result_labels(qual),{' '},'(',result_units(qual),')');
        numberOfXTicks = 4;
        xData = get(h,'XData');
        set(gca,'Xtick',logspace(log(xData(2))/log(10),log(xData(end))/log(10),numberOfXTicks))
        if qual==9
            caxis([0 3500])
        elseif qual==8
            caxis([0 3500])
        elseif qual==7
            caxis([0 0.85])
        end
        
    end
%     h = text(-0.25, 0.5, 'row 2');
    saveas(gcf,strcat(savedir,'Dc_Db_K_b_t','_',save_labels{qual},...
        '_hm.eps'),'epsc')
end
%% Exploring the effects of a stationary buffer
% Load PS1
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=cell(size(un{4}));
resrange=find(results(:,5)==un{5}(6));
resrange=intersect(find(results(:,6)==un{6}(6)),resrange);
resrange=intersect(find(results(:,2)==un{2}(1)),resrange);
resrangepreloop=(resrange);
numd=3;
numa=4;
for qual=7:9
    figure('units','normalized','outerposition',[0 0 1 1])
    for ii=1:numpervar
        subplot(numd,numa,ii)
        hold on
        resrangel=intersect(find(results(:,1)==un{1}(ii)),resrangepreloop);
        for jj=1:numpervar
            resrangeK{jj}=intersect(find(results(:,4)==un{4}(jj)),resrangel);
            legText{jj}=strcat('K=',num2str(round(un{4}(jj),2,'significant')),'\muM');
            x = results(resrangeK{jj},3);
            dc= results(resrangeK{jj},qual);
            if jj<8
                h=plot(x,dc,'-o','MarkerSize',5,'LineWidth',2);
            else
                h=plot(x,dc,'-x','MarkerSize',5,'LineWidth',2);
            end
        end
        set(gca,'XScale','log')
        title(strcat('D_c=',num2str(round(un{1}(ii),1,'significant'))))
%         if qual==9
%             axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
%             300 2200])
%         elseif qual==8
%             axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
%             0 320])
%         elseif qual==7
%             axis([min(results(resrangepreloop,3)) max(results(resrangepreloop,3)) ...
%             0 1])
%         end
        ylabel(strcat(result_labels(qual),{' '},'(',result_units(qual),')'))%,'FontSize',20)
        xlabel(strcat(result_labels(3),{' '},'(',result_units(3),')'))%,'FontSize',20)
        

%         numberOfXTicks = 3;
%         xData = get(h,'XData');
%         set(gca,'Xtick',logspace(log(xData(1))/log(10),log(xData(end))/log(10),numberOfXTicks))
    end
    if qual==7
        legend(legText,'Location','best','Box','off','FontSize',12)
    else
        legend(legText,'Location','best','Box','off','FontSize',12)
    end
%         saveas(gcf,strcat(savedir,'ch3_sbuffers-',num2str(pos),save_labels{qual},...
%         '.eps'),'epsc')
end
%% Chapter fig for stationary buffers (they're pretty similar in the expected range of Dc)
% Load PS1
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
title_labels={'D_c','D_b','b_t','K','NPC','a','Ca^{2+} amplitude','TTP','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);

for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=title_labels;
% resrange=find(results(:,5)==un{5}(6));
% resrange=intersect(find(results(:,6)==un{6}(6)),resrange);
% resrange=intersect(find(results(:,2)==un{2}(1)),resrange);
% resrange=intersect(find(results(:,1)==un{1}(5)),resrange);
% resrangepreloop=(resrange);
numpervar=sqrt(size(results,1));
resrangepreloop=1:numpervar^2;
numd=3;
numa=4;
baseSize=20;
cl=[63,0,125]/255;
ch = [252,251,253]/255;
for qual=7:9
    figure%('units','normalized','outerposition',[0 0 1 1])
        hold on
        pmap=((1:64)'*ch+(64-(1:64)')*cl)/64;
        colormap(pmap)
        count=0;
        for jj=1:round(numpervar/7):numpervar
            count=count+1;
            resrangeK{count}=intersect(find(results(:,4)==un{4}(jj)),resrangepreloop);
            legText{count}=strcat('K=',num2str(round(un{4}(jj),2,'significant')),'\muM');
            x = results(resrangeK{count},3);
            dc= results(resrangeK{count},qual);
            plot(x,dc,'LineWidth',3,'Color',(jj*ch+(numpervar-jj)*cl)/numpervar);
        end
        set(gca,'XScale','log')
        axis('tight')
        set(gca,'FontSize',baseSize*0.8)
        title(strcat('Effect of a stationary buffer on',{' '},title_labels{qual}),'FontSize',1.2*baseSize)
        ylabel(strcat(result_labels(qual),{' '},'(',result_units(qual),')'),'FontSize',baseSize)
        xlabel(strcat(result_labels(3),{' '},'(',result_units(3),')'),'FontSize',baseSize)
%         colorbar('Yscale', 'log')
        c=colorbar;
        set(gca,'ColorScale','log')
        caxis([1e-6,max(un{4})])
        c.Label.String = strcat(result_labels(4),{' '},'(',result_units(4),')');
        c.Label.FontSize = baseSize;
%     if qual==7
%         legend(legText,'Location','best','Box','off','FontSize',0.8*baseSize)
%     else
%         legend(legText,'Location','best','Box','off','FontSize',0.8*baseSize)
%     end
%         saveas(gcf,strcat(savedir,'ch3_sbuffer-',num2str(pos),'_D_c=',num2str(round(un{1}(5),2,'significant')),'_',save_labels{qual},...
%         '.eps'),'epsc')
    saveas(gcf,strcat(savedir,'ch3_sbuffer_sh-',num2str(pos),'_D_c=00176_',save_labels{qual},...
        '.eps'),'epsc')
end

%% Chapter fig (heatmap) for stationary buffers (they're pretty similar in the expected range of Dc)
% New! Get ps from huia 
% Load parPS_bK12.mat
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
title_labels={'D_c','D_b','b_t','K','NPC','a','Ca^{2+} amplitude','TTP','FDHM','ttbase'};
% numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);

for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
legText=title_labels;
% resrange=find(results(:,5)==un{5}(6));
% resrange=intersect(find(results(:,6)==un{6}(6)),resrange);
% resrange=intersect(find(results(:,2)==un{2}(1)),resrange);
% resrange=intersect(find(results(:,1)==un{1}(5)),resrange);
% resrangepreloop=(resrange);
numpervar=sqrt(size(results,1));
resrangepreloop=1:numpervar^2;

numd=3;
numa=4;
baseSize=42;
for qual=7:9
    fpos=[400 250 600 500];
    figure('pos',fpos)%('units','normalized','outerposition',[0 0 1 1])
    ax=gca;
    ax.Position=[0.21 0.25 0.6 0.6];hold on
    x = unique(results(resrangepreloop,3));
    y = unique(results(resrangepreloop,4));
    dc= reshape(results(resrangepreloop,qual),[numpervar,numpervar]);
    surf(x,y,dc,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    view(0,90)
%     axis('square')
    axis('tight')
    set(gca,'FontSize',baseSize*0.8)
%     title(strcat('Effect of a stationary buffer on',{' '},title_labels{qual}),'FontSize',baseSize,'FontWeight','normal')
    ylabel(strcat(result_labels(4),{' '},'(',result_units(4),')'),'FontSize',baseSize*0.9)
    xlabel(strcat(result_labels(3),{' '},'(',result_units(3),')'),'FontSize',baseSize*0.9)
	c=colorbar;
    c.Label.String = strcat(result_labels(qual),{' '},'(',result_units(qual),')');
    c.Label.FontSize = baseSize*0.8;
    cpos=c.Label.Position;
    c.Label.Position=[cpos(1)-1e-2 cpos(2:end)];
    axlims=get(gcf,'CurrentAxes');
    newx=10.^(floor(log10(axlims.XLim(1))):floor(log10(axlims.XLim(end))));
    set(gca,'XTick',newx(1:round(end/3):end));
%     if qual==7
%         legend(legText,'Location','best','Box','off','FontSize',0.8*baseSize)
%     else
%         legend(legText,'Location','best','Box','off','FontSize',0.8*baseSize)
%     end
        saveas(gcf,strcat(savedir,'ch3_sbuffer_hm-',num2str(pos),'_D_c=00174','_',save_labels{qual},...
        '.eps'),'epsc')
end

%% Compare different transient lengths
c1=[27 158 119]/255;
c2=[215 94 0]/255;
c3=[117 100 175]/255;


[cutca,~,~]=nonlinearDiff(1e-2,[],fits,crtpts,2.6,3000);
[cutcall,~,~]       =nonlinearDiff_longer(1e-2,[],0.8,fits,crtpts,2.6,3000);
[cutcals,~,t]       =nonlinearDiff_longer(1e-2,[],1.2,fits,crtpts,2.6,3000);
figure
hold on
h1=plot(t,cutcall(:,end),'--','LineWidth',2,'Color',c2);
h3=plot(t,cutcall(:,1),'LineWidth',2,'Color',c2);
h5=plot(t,cutcals(:,end),'--','LineWidth',2,'Color',c1);
h6=plot(t,cutcals(:,1),'LineWidth',2,'Color',c1);
h2=plot(t,cutca(:,end),'--','LineWidth',2,'Color',c3);
h4=plot(t,cutca(:,1),'LineWidth',2,'Color',c3);
set(gca,'FontSize',0.7*baseSize)
legend([h2,h4,h1,h3,h5,h6],{'normal cyt','normal nuc','longer cyt','longer nuc','shorter cyt','shorter nuc'},...
    'FontSize',0.7*baseSize,'Box','off')
xlabel('Time (ms)','FontSize',0.8*baseSize)
ylabel('[Ca^{2+}]','FontSize',0.8*baseSize)
% title('D_c=0.01\mum^2/ms','FontSize',0.9*baseSize,'FontWeight','normal')
axis([0 3000 0 1])
saveas(gcf,strcat(savedir,'longer_trans_comparison.eps'),'epsc')

%% Compare time to reach half max with different transient lengths
[cutcaln,~,~]   =nonlinearDiff_longer(1e-2,[],1,fits,crtpts,2.6,3000);
[cutcall,~,~]   =nonlinearDiff_longer(1e-2,[],0.8,fits,crtpts,2.6,3000);
[cutcals,r,t]   =nonlinearDiff_longer(1e-2,[],1.2,fits,crtpts,2.6,3000);
hmax=(max(cutcaln(:,end))+min(cutcaln(:,end)))/2;
thmax=zeros(3,size(cutcaln,2));
for ii=1:size(cutcaln,2)
   thmax(1,ii)= range(find(cutcaln(:,ii)>=hmax));
   thmax(2,ii)= range(find(cutcall(:,ii)>=hmax));
   thmax(3,ii)= range(find(cutcals(:,ii)>=hmax));
end
figure
plot(r,thmax./thmax(1,:))
%% Compare different transient heights
[cutca,~,~]=nonlinearDiff(1e-2,[],fits,crtpts,2.6,3000);
[cutcahs,~,~]       =nonlinearDiff_higher(1e-2,[],0.8,fits,crtpts,2.6,3000);
[cutcahh,~,t]       =nonlinearDiff_higher(1e-2,[],1.2,fits,crtpts,2.6,3000);
figure
hold on
h1=plot(t,cutcahs(:,end),'--','LineWidth',2,'Color',c2);
h3=plot(t,cutcahs(:,1),'LineWidth',2,'Color',c2);
h5=plot(t,cutcahh(:,end),'--','LineWidth',2,'Color',c1);
h6=plot(t,cutcahh(:,1),'LineWidth',2,'Color',c1);
h2=plot(t,cutca(:,end),'--','LineWidth',2,'Color',c3);
h4=plot(t,cutca(:,1),'LineWidth',2,'Color',c3);
legend([h2,h4,h1,h3,h5,h6],{'normal cyt','normal nuc','shorter cyt','shorter nuc','teller cyt','taller nuc'},'FontSize',14)
xlabel('Time (ms)','FontSize',18)
ylabel('[Ca^{2+}]','FontSize',18)
title('D_c=0.01\mum^2/ms','FontSize',20)
set(gca,'FontSize',12)
axis([0 3000 0 1.2])

%% Chapter Fig 13?
% Add lt and ht
numEachVar=11;
numSteps = numEachVar^6;
ns = numEachVar-1;
param_type=1;
ht=0.5+1*repmat(reshape(repmat((0:(1/ns):1),[numEachVar^2,1]),numEachVar^3,1),[numEachVar^3,1]);
lt=0.5+1*repmat(reshape(repmat((0:(1/ns):1),[numEachVar,1]),numEachVar^2,1),[numEachVar^4,1]);
results=[results,ht,lt];
results(:,12)=1./results(:,12);
%% Effect of changes to max and Fdhm of cyt input on quals - heatmap
%load('/Users/hhunt1/Documents/huia/parPS_hl1')
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc} (\muM)','TTP (nuc) (ms)','FDHM (nuc) (ms)','Time to base','Max [Ca^{2+}]_{cyt}','FDHM (cyt)'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms','\muM','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase','hcyt','lcyt'};
numpervar=11;%round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=[1:6,11:12]
    un{unvar}=unique(results(:,unvar));
end
resrange=find(results(:,1)==un{1}(4));
resrangepreloop=intersect(find(results(:,2)==un{2}(1)),resrange); % Pref 2
resrangeb=intersect(find(results(:,2)==un{2}(5)),resrange);
resrangeb=intersect(find(results(:,3)==un{3}(2)),resrangeb);
resrangeb=intersect(find(results(:,4)==un{4}(6)),resrangeb);
for fixvar=[3,4]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(1)),resrangepreloop);
end
for qual=7:9
    ii=12;
    jj=11;
    
    figure('position',[0,200,560,450])
    if 1
    dc=reshape(results(resrangepreloop,qual),[numpervar,numpervar]);
    x=unique(results(resrangepreloop,ii),'stable')*291.2;
    y=unique(results(resrangepreloop,jj),'stable');
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
%     set(gca,'XScale','log')
    else
    x = results(resrangepreloop,ii);
    y = results(resrangepreloop,jj);
    dc= results(resrangepreloop,qual);
    scatter(x,y,4500,dc,'filled')
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    
%     axis([min(x-2) max(x+2) min(y-0.25) max(y+0.25)])
%     set(gca,'XScale','log')
    end
    c=colorbar;
    c.Label.String=strcat(result_labels{qual});
    c.Label.FontSize=0.7*sqSize;
    ax=gca;
    ax.Position=[0.2 0.25 0.6 0.6];
    set(gca,'FontSize',sqSize*0.8)
%     title(strcat(result_labels{qual}),'FontSize',sqSize*0.7,'FontWeight','normal')
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',sqSize*0.8)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',sqSize*0.8)
    saveas(gcf,strcat(savedir,'chapterfig_cyttrans-',num2str(pos),'_',save_labels{qual},'_','nobuf','.eps'),'epsc')
end
%% Effect of changes to max and Fdhm of cyt input on quals - heatmap % change
%load('/Users/hhunt1/Documents/huia/parPS_hl1')
savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base','Max [Ca^{2+}]_{cyt}','FDHM (cyt)'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms','\muM','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase','hcyt','lcyt'};
numpervar=11;%round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=[1:6,11:12]
    un{unvar}=unique(results(:,unvar));
end
resrange=find(results(:,1)==un{1}(4));
resrangepreloop=intersect(find(results(:,2)==un{2}(1)),resrange); % Pref 2
resrangeb=intersect(find(results(:,2)==un{2}(5)),resrange);
resrangeb=intersect(find(results(:,3)==un{3}(2)),resrangeb);
resrangeb=intersect(find(results(:,4)==un{4}(6)),resrangeb);
for fixvar=[3,4]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(1)),resrangepreloop);
end
% Find result for no change in Max or FDHM
fixresrange=intersect(find(results(:,11)==1),resrangepreloop);
fixresrange=intersect(find(results(:,12)==1),fixresrange);
for qual=7:9
    ii=12;
    jj=11;
    
    figure('position',[0,200,560,450])
    if 1
    dc=100*reshape(results(resrangepreloop,qual),[numpervar,numpervar])/results(fixresrange,qual)-100;
    x=unique(results(resrangepreloop,ii),'stable')*100-100;
    y=unique(results(resrangepreloop,jj),'stable')*100-100;
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
%     set(gca,'XScale','log')
    else
    x = results(resrangepreloop,ii)*100-100;
    y = results(resrangepreloop,jj)*100-100;
    dc= results(resrangepreloop,qual);
    scatter(x,y,4500,dc,'filled')
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
%     axis([min(x-2) max(x+2) min(y-0.25) max(y+0.25)])
%     set(gca,'XScale','log')
    end
    title(strcat(result_labels{qual},{' '},'(% increase)'),'FontSize',25)
    ylabel(strcat(result_labels(jj),{' '},'(% increase)'),'FontSize',20)
    xlabel(strcat(result_labels(ii),{' '},'(% increase)'),'FontSize',20)
    axis([-20 20 -20 20 ])
    if qual==9
        caxis([-15,15])
    elseif qual==8
        caxis([-3 3])
    else
        caxis([-20 20])
    end
    colorbar
    if saveBool
    saveas(gcf,strcat(savedir,'chapterfig_cyttrans_pc-',num2str(pos),'_',save_labels{qual},'_','nobuf','.eps'),'epsc')
    end
end

%% buffers
for qual=7:9
    ii=12;
    jj=11;
    
    figure('position',[0,200,560,450])
    if 1
    dc=reshape(results(resrangeb,qual),[numpervar,numpervar]);
    x=unique(results(resrangeb,ii),'stable');
    y=unique(results(resrangeb,jj),'stable');
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    view(90,90)
    axis([min(x) max(x) min(y) max(y)])
%     set(gca,'XScale','log')
    else
    x = results(resrangeb,ii);
    y = results(resrangeb,jj);
    dc= results(resrangeb,qual);
    scatter(x,y,4500,dc,'filled')
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    axis([min(x-2) max(x+2) min(y-0.25) max(y+0.25)])
%     set(gca,'XScale','log')
    end
    title(strcat(result_labels{qual}),'FontSize',1.25*baseSize)
    ylabel(strcat(result_labels(jj)),'FontSize',baseSize)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize)
    colorbar
    saveas(gcf,strcat(savedir,'chapterfig_cyttrans-',num2str(pos),'_',save_labels{qual},'_','cam','.eps'),'epsc')
end
%% Chapter Fig 14?
% Effect of cyt fdhm vs qual. With diff lines for diff Dc
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base','Max [Ca^{2+}]_{cyt}','FDHM (cyt)'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms','\muM','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase','hcyt','lcyt'};
numpervar=11;%round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=[1:6,11:12]
    un{unvar}=unique(results(:,unvar));
end
resrange=1:size(results,1);
resrangepreloop=intersect(find(results(:,2)==un{2}(1)),resrange); % Pref 2
resrangeb=intersect(find(results(:,2)==un{2}(5)),resrange);
resrangeb=intersect(find(results(:,3)==un{3}(2)),resrangeb);
resrangeb=intersect(find(results(:,4)==un{4}(6)),resrangeb);
for fixvar=[3,4]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(1)),resrangepreloop);
end
for lvar=[11,12]    
resrangeloop=intersect(find(results(:,lvar)==un{lvar}(6)),resrangepreloop);
figure
numa=2;numd=2;
for qual=7:9
    subplot(numd,numa,qual-6)
    hold on
    for ii=1:size(un{1},1)
       resrangel=intersect(find(results(:,1)==un{1}(ii)),resrangeloop);
       legText{ii}=strcat('D_c=',num2str(round(un{1}(ii),1,'significant')),'\mum^2/ms');
       x=results(resrangel,23-lvar);
       dc=results(resrangel,qual);
       if ii<7
           plot(x,dc,'LineWidth',2)
       else
%            plot(x,dc,'--','LineWidth',2)
       end
    end
    set(gca,'FontSize',16)
%     title(strcat(result_labels{qual}),'FontSize',25)
    xlabel(strcat(result_labels(23-lvar),'(',result_units(23-lvar),')'),'FontSize',18)
    ylabel(strcat(result_labels(qual),{' '},'(',result_units(qual),')'),'FontSize',18)
end
    legend(legText,'Location','best','box','off')
end

%% Chapter Fig 15?
% Effect of cyt fdhm vs qual. Heatmaps
% load('/Users/hhunt1/Documents/huia/parPS_hl1')
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base','Max [Ca^{2+}]_{cyt}','FDHM_{cyt}'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms','\muM','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase','hcyt','lcyt'};
numpervar=11;%round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=[1:6,11:12]
    un{unvar}=unique(results(:,unvar));
end
resrange=1:size(results,1);
resrangepreloop=intersect(find(results(:,2)==un{2}(1)),resrange); % Pref 2
resrangeb=intersect(find(results(:,2)==un{2}(5)),resrange);
resrangeb=intersect(find(results(:,3)==un{3}(2)),resrangeb);
resrangeb=intersect(find(results(:,4)==un{4}(6)),resrangeb);
for fixvar=[3,4]
    resrangepreloop=intersect(find(results(:,fixvar)==un{fixvar}(1)),resrangepreloop);
end

for lvar=[11,12]    
resrangeloop=intersect(find(results(:,lvar)==un{lvar}(6)),resrangepreloop);
% figure
numa=2;numd=2;
for qual=7:9
    figure
    y=unique(results(resrangeloop,1));
    x=unique(results(resrangeloop,23-lvar));
    dc=reshape(results(resrangeloop,qual),[numpervar,numpervar]);
    surf(x,y,dc,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(0,90)
    ax=gca;
    ax.Position=[0.22 0.25 0.6 0.6];
    set(gca,'YScale','log')
    set(gca,'FontSize',sqSize*0.8)
    axlims=get(gcf,'CurrentAxes');
    newy=10.^(floor(log10(axlims.YLim(1))):floor(log10(axlims.YLim(end))));
    set(gca,'YTick',newy(1:round(end/3):end));
    ytickformat('%.0e');
	axis('tight')
    xl=xlabel(strcat(result_labels(23-lvar),{' '},'(',result_units(23-lvar),')'),'FontSize',sqSize*0.8);
    yl=ylabel(strcat(result_labels(1),{' '},'(',result_units(1),')'),'FontSize',sqSize*0.8);
    ypos=get(yl,'Position');
    set(yl,'Position',[ypos(1)+1e-2 ypos(2:end)])
    xpos=get(xl,'Position');
    set(xl,'Position',[xpos(1) xpos(2)+1e-5 xpos(3:end)])
    c=colorbar;
    c.Label.String=strcat(result_labels(qual),{' '},'(',result_units(qual),')');
    c.Label.FontSize=sqSize*0.8;
    saveas(gcf,string(strcat(savedir,'lh_buffers-',save_labels(23-lvar),'_',save_labels(qual),'.eps')),'epsc')
end
end
%% Id buffer
%% Investigate individual buffers
% d_c v D_b v fdhm 
% Load PS1
savedir='/Volumes/UniWD/AnaNucleus/numPS/';
savebool=0;
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
    savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
end
pos=1;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'\mum^2/ms','\mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','Max','ttpeak','FDHM','ttbase'};
numpervar=round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=1:6
    un{unvar}=unique(results(:,unvar));
end
resrange=find(results(:,5)==un{5}(5));
resrangepreloop=intersect(find(results(:,6)==un{6}(5)),resrange);
% fixval=[1,1,11,8,1,1];
buffer_label={'Calmodulin','NE','Calreticulin Site C','Calreticulin Site P','ATP'};
fixval=cell(1,5);
fixval{1}=[1,1,1,6,1,1];
fixval{2} =[1,1,5,4,1,1];
fixval{3}=[1,1,9,9,1,1];
fixval{4}=[1,1,7,8,1,1];
fixval{5}=[1,1,8,9,1,1];
bufD=[2,0,0.027,0.027,0.14];

minfdhm=min(dataStruct(cellit).nfdhm)*0.8;
maxfdhm=max(dataStruct(cellit).nfdhm)*1.2;
minttp=min(dataStruct(cellit).nttpeak)*0.8;
maxttp=max(dataStruct(cellit).nttpeak)*1.2;

for bufit=1%1:5%[1,3:5]
    resrangeloop=resrangepreloop;
for fixvar=[3,4]
    resrangeloop=intersect(find(results(:,fixvar)==un{fixvar}(fixval{bufit}(fixvar))),resrangeloop);
end
if bufit==2
    resrangeloop=intersect(find(results(:,2)==un{2}(1)),resrangeloop);
end

if bufit==2
    dcfdhm=results(resrangeloop,9);
    dcttp=results(resrangeloop,8);
    sweetspot=(dcfdhm>minfdhm&dcfdhm<maxfdhm)+(dcttp>minttp&dcttp<maxttp);
else
dcfdhm=reshape(results(resrangeloop,9),[numpervar,numpervar]);
dcttp=reshape(results(resrangeloop,8),[numpervar,numpervar]);
sweetspot=2*(dcfdhm>minfdhm&dcfdhm<maxfdhm)+(dcttp>minttp&dcttp<maxttp);
end
% mesh(sweetspot)
ii=1;
jj=2;

figure('position',[0,200,560,450])
x=unique(results(resrangeloop,ii));
y=unique(results(resrangeloop,jj));
xn=results(resrangeloop,ii);
yn=results(resrangeloop,jj);
if bufit==2
    plot(x,sweetspot,'LineWidth',2)
else
surf(x,y,sweetspot','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
view(90,90)
axis([min(x) max(x) min(y) max(y)])

vertex_matrix = [min(x) bufD(bufit) 1; ...
                    min(x) bufD(bufit) 3; ...
                    max(x) bufD(bufit) 1; ...
                    max(x) bufD(bufit) 3];
    %             faces_matrix = [4*ii-3 4*ii-2 4*ii 4*ii-1];
faces_matrix = [1 2 4 3];
patch('Vertices',vertex_matrix,'Faces',faces_matrix, ...
    'FaceColor','y','EdgeColor', ...
    c2,'FaceAlpha',0.2,'LineWidth',3)
if bufit==1
    vertex_matrix = [min(x) 40 1; ...
                    min(x) 40 3; ...
                    max(x) 40 1; ...
                    max(x) 40 3];
    %             faces_matrix = [4*ii-3 4*ii-2 4*ii 4*ii-1];
faces_matrix = [1 2 4 3];
patch('Vertices',vertex_matrix,'Faces',faces_matrix, ...
    'FaceColor','y','EdgeColor', ...
    c2,'FaceAlpha',0.2,'LineWidth',3)
end

axis('square')
set(gca,'YScale','log')
% clabel(C,h,'manual','Color','w','FontSize',baseSize*0.8)
% colorbar
ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',baseSize)
end
set(gca,'XScale','log')
set(gca,'FontSize',baseSize*0.8)
title(strcat(buffer_label{bufit}),'FontSize',1.25*baseSize)
xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',baseSize)
% saveas(gcf,strcat(savedir,'chapterfig_buffers-',num2str(bufit),'.eps'),'epsc')
end

%% There is no sweetspot
tol = 0;
minfdhm=min(dataStruct(cellit).nfdhm)*(1-tol);
maxfdhm=max(dataStruct(cellit).nfdhm)*(1+tol);
minttp=min(dataStruct(cellit).nttpeak)*(1-tol);
maxttp=max(dataStruct(cellit).nttpeak)*(1+tol);

% dcfdhm=results(:,9);
% dcttp=results(:,8);
dcfdhm=results(:,10);
dcttp=results(:,9);
sweetspot=(dcfdhm>minfdhm&dcfdhm<maxfdhm)+(dcttp>minttp&dcttp<maxttp);
closeres=find(sweetspot);
possSpots=sum(find(sweetspot>1))
%
figure
hold on
% plot(dcttp,dcfdhm,results(:,7),'x','MarkerSize',10,'LineWidth',2)
% scatter(results(:,8),results(:,9),10,results(:,7),'LineWidth',2)
scatter(results(:,9),results(:,10),10,results(:,8),'LineWidth',2)
plot([minttp minttp],[minfdhm maxfdhm],'r')
plot([maxttp maxttp],[minfdhm maxfdhm],'r')
plot([minttp maxttp],[maxfdhm maxfdhm],'r')
plot([minttp maxttp],[minfdhm minfdhm],'r')
% plot([minttp minttp],[minfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([maxttp maxttp],[minfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([minttp maxttp],[maxfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([minttp maxttp],[minfdhm minfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% axis([minttp*0.7 maxttp*1.4 minfdhm*0.7 maxfdhm*1.3])

%% There is no sweetspot
%load('cell7-1')
%load(parpshl - reshl
nfdhm=dataStruct(cellit).nfdhm;
fdhmerr=std(nfdhm);
meanfdhm=mean(nfdhm);

nttp=dataStruct(cellit).nttpeak;
ttperr=std(nttp);
meanttp=mean(nttp);

% Find optimal frontiers
% parPS1, change ttp
%restrict ot results with fdhm between 0 and 1500 and ttp between 0 and 700
if 0
resres=results(results(:,9)<1500,:);
resres=resres(resres(:,8)<700,:);
resres=resres(resres(:,1)<1,:);
resres=resres(resres(:,2)<1,:);
resres=resres(resres(:,5)==1.05,:);
resres=resres(resres(:,6)==1.96,:);
% ttp
% [~,fdhmsort]=sort(results(:,9));
topttpx=[];
topttpy=[];
botttpx=[];
botttpy=[];
for ii=[0:100:1400]
    temprr=resres(resres(:,9)>ii,:);
    [tempval,temploc]=min(temprr(:,8));
    topttpx=[topttpx,tempval];
    topttpy=[topttpy,temprr(temploc,9)];
    
    temprr=resres(resres(:,9)<ii&resres(:,9)>0,:);
    if ~isempty(temprr)
    [tempval,temploc]=max(temprr(:,8));
    botttpx=[botttpx,tempval];
    botttpy=[botttpy,temprr(temploc,9)];
    end
end
end
% parpshl
% load(parpshl), change ttp
resreshl=reshl(reshl(:,9)<1500,:);
resreshl=resreshl(resreshl(:,8)<700,:);
resreshl=resreshl(resreshl(:,1)<1,:);
resreshl=resreshl(resreshl(:,2)<1,:);
resreshl=resreshl(resreshl(:,12)==2,:);
% ttp
% [~,fdhmsort]=sort(results(:,9));
topttpxhl=[];
topttpyhl=[];
botttpxhl=[];
botttpyhl=[];
for ii=[0:100:1400]
    temprr=resreshl(resreshl(:,9)>ii,:);
    [tempval,temploc]=min(temprr(:,8));
    topttpxhl=[topttpxhl,tempval];
    topttpyhl=[topttpyhl,temprr(temploc,9)];
    
    temprr=resreshl(resreshl(:,9)<ii&resreshl(:,9)>0,:);
    if ~isempty(temprr)
    [tempval,temploc]=max(temprr(:,8));
    botttpxhl=[botttpxhl,tempval];
    botttpyhl=[botttpyhl,temprr(temploc,9)];
    end
end

baseSize=38;
fpos=[400 250 760 340];
figure('pos',fpos)
hold on
l1=patch([topttpx,botttpx(end:-1:1)],[topttpy,botttpy(end:-1:1)],[0 0 1],'FaceAlpha',0.2,'EdgeColor','none');
l2=patch([topttpxhl,botttpxhl(end:-1:1)],[topttpyhl,botttpyhl(end:-1:1)],[1 0 0],'FaceAlpha',0.2,'EdgeColor','none');
errorbar(meanttp,meanfdhm,fdhmerr,'or','LineWidth',2,'MarkerSize',10,'CapSize',10);
l3=errorbar(meanttp,meanfdhm,ttperr,'horizontal','or','LineWidth',2,'MarkerSize',10,'CapSize',10);
legend([l1,l2,l3],{'original FDHM','increased FDHM','cell data'},'Location','northwest')
set(gca,'FontSize',baseSize*0.7)
xlabel('TTP (ms)','FontSize',0.8*baseSize)
ylabel('FDHM (ms)','FontSize',0.8*baseSize)
axis([0 700 0 1500])
%% Closest to a sweetspot
tol = 0.2;
minfdhm=min(dataStruct(cellit).nfdhm)*(1-tol);
maxfdhm=max(dataStruct(cellit).nfdhm)*(1+tol);
minttp=min(dataStruct(cellit).nttpeak)*(1-tol);
maxttp=max(dataStruct(cellit).nttpeak)*(1+tol);

dcfdhm=results(:,9);
dcttp=results(:,8);
sweetspot=2*(dcfdhm>minfdhm&dcfdhm<maxfdhm)+3*(dcttp>minttp&dcttp<maxttp);
find(sweetspot==5);
%% Surf simulation
% load('cell7-1')
savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
% cellit=6;
cellit=1;
baseSize=42;
D=0.0176;
timeInt=10;%(30*dataStruct(cellit).timeStep);
timeEnd=(1e4*dataStruct(cellit).timeStep);
time=0:timeInt:timeEnd;
fits=dataStruct(cellit).fstfit{1};
crtpts=dataStruct(cellit).crtpts{1};
[numSol,~,~]=nonlinearDiff(D,[],fits,crtpts,timeInt,timeEnd);
timevec=0:timeInt:timeEnd;
spacevec=linspace(0,2,2e2+1);
tpart=19:100;
figure('position',[0,200,600,500])
surf([spacevec,2+spacevec(2:end)],timevec(tpart)-timevec(tpart(1)),numSol(tpart,[end:-1:1 2:end]),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
view(90,90)
set(gca,'FontSize',baseSize*0.7)
ylabel('time (ms)','FontSize',baseSize*0.8)
xlabel('position (\mum)','FontSize',baseSize*0.8)
ax=gca;
ax.Position=[0.22 0.25 0.6 0.6];
c=colorbar;
c.Label.String='[Ca^{2+}] (\muM)';
c.Label.FontSize=baseSize*0.8;
axis('tight')
saveas(gcf,strcat(savedir,'sim_1D_D-00176_cell-',num2str(cellit),'.eps'),'epsc')
if 1
figure
surf([spacevec(1:10:end),2+spacevec(1:10:end)],timevec(tpart)-timevec(tpart(1)),numSol(tpart,[end:-10:1 1:10:end]),'FaceColor','interp','FaceLighting','gouraud');
view(45,50)
set(gca,'FontSize',baseSize*0.7)
ylabel('time (ms)','FontSize',baseSize*0.8)
xl=xlabel('position (\mum)','FontSize',baseSize*0.8);
c=colorbar;
c.Label.String='[Ca^{2+}] (\muM)';
c.Label.FontSize=baseSize*0.8;
axis('tight')
ax=gca;
ax.Position=[0.1 0.25 0.6 0.6];
xpos=get(xl,'Position');
set(xl,'Position',[xpos(1)+2 xpos(2)+2e2 xpos(3:end)])
saveas(gcf,strcat(savedir,'sim_mesh_1D_D-00176_cell-',num2str(cellit),'.eps'),'epsc')
end

%% Find cell that will generate fdhm and ttp of data
if 0
    load('/Users/hhunt1/Documents/huia/parPS_hl1')
    load('cell7-1')
end
% resrange=find(results(:,5)==1.05);
% resrange2=intersect(find(results(:,6)==1.96),resrange);
resrange2=1:size(results,1);
resPlot=results(resrange2,:);
resPlot(:,8)=resPlot(:,8)*10-182;
resPlot(resPlot(:,8)<0,:)=[];
resplot2=resPlot(resPlot(:,12)==2,:);
tol = 0;
minfdhm=min(dataStruct(cellit).nfdhm)*(1-tol);
maxfdhm=max(dataStruct(cellit).nfdhm)*(1+tol);
minttp=min(dataStruct(cellit).nttpeak)*(1-tol);
maxttp=max(dataStruct(cellit).nttpeak)*(1+tol);

fpos=[400 250 760 340];
figure('pos',fpos)
hold on
plot(resPlot(:,8),resPlot(:,9),'x','MarkerSize',10,'LineWidth',2)
plot([minttp minttp],[minfdhm maxfdhm],'r','LineWidth',2)
plot([maxttp maxttp],[minfdhm maxfdhm],'r','LineWidth',2)
plot([minttp maxttp],[maxfdhm maxfdhm],'r','LineWidth',2)
plot([minttp maxttp],[minfdhm minfdhm],'r','LineWidth',2)
% plot([minttp minttp],[minfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([maxttp maxttp],[minfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([minttp maxttp],[maxfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([minttp maxttp],[minfdhm minfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% axis([minttp*0.7 maxttp*1.4 minfdhm*0.7 maxfdhm*1.3])

set(gca,'FontSize',baseSize*0.7)
xlabel('TTP (ms)','FontSize',0.8*baseSize)
ylabel('FDHM (ms)','FontSize',0.8*baseSize)
axis([0 1000 0 2000])
% saveas(gcf,strcat(savedir,'hl_buf_fdhm_ttp.eps'),'epsc')
%% Plot solution to match data
% load('cell7-1')
timeInt=10;%(30*dataStruct(cellit).timeStep);
timeEnd=(1e4*dataStruct(cellit).timeStep);
time=0:timeInt:timeEnd;
% [numSol]=nonlinearDiff_hl(0.0251188643150958,[0.00630957344480194,100,10,1,2],1,2,fits,crtpts,timeInt,timeEnd);
[numSol]=nonlinearDiff_hl(0.000398107170553497,[0.1,1,10,1,2],1,2,fits,crtpts,timeInt,timeEnd);

%% long sim buf hl
sqSize=48;
fpos=[400 250 760 340];
figure('pos',fpos)
surf((0:198)*4/198,time-182,numSol(:,[end:-1:1 2:end]),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
axis('tight')
axis([0 4 0 1000])
view(90,90)
c=colorbar;
set(gca,'FontSize',sqSize*0.7)
ylabel('time (ms)','FontSize',sqSize*0.8)
xlabel('position (\mum)','FontSize',sqSize*0.8)
cpos=c.Position;
c.Position=[cpos(1) cpos(2:end)];
c.FontSize=0.8*sqSize;
set(get(c,'Title'),'String','      (\muM)','FontSize',sqSize*0.7)
% c.Label.String='[Ca^{2+}] (\muM)';
% c.Label.FontSize=0.7*baseSize;
clpos=c.Title.Position;
c.Title.Position=[clpos(1)+20 clpos(2)-30 clpos(3:end)];
caxis([0 0.7])
saveas(gcf,strcat(savedir,'match_data_sim.eps'),'epsc')

%% long sim just diff
cellit=1;
D=0.0176;
timeInt=10;%(30*dataStruct(cellit).timeStep);
timeEnd=(1e4*dataStruct(cellit).timeStep);
time=0:timeInt:timeEnd;
fits=dataStruct(cellit).fstfit{1};
crtpts=dataStruct(cellit).crtpts{1};
[numSol,~,~]=nonlinearDiff(D,[],fits,crtpts,timeInt,timeEnd);

sqSize=48;
fpos=[400 250 760 340];
figure('pos',fpos)
surf((0:400)*4/400,time-182,numSol(:,[end:-1:1 2:end]),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
axis('tight')
axis([0 4 0 1000])
view(90,90)
c=colorbar;
set(gca,'FontSize',sqSize*0.7)
ylabel('time (ms)','FontSize',sqSize*0.8)
xlabel('position (\mum)','FontSize',sqSize*0.8)
cpos=c.Position;
c.Position=[cpos(1) cpos(2:end)];
c.FontSize=0.8*sqSize;
set(get(c,'Title'),'String','      (\muM)','FontSize',sqSize*0.7)
% c.Label.String='[Ca^{2+}] (\muM)';
% c.Label.FontSize=0.7*baseSize;
clpos=c.Title.Position;
c.Title.Position=[clpos(1)+20 clpos(2)-270 clpos(3:end)];
saveas(gcf,strcat(savedir,'match_data_sim_diff.eps'),'epsc')