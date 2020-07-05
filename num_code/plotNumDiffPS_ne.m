%% Plot numDiffPS_ne results

% Gridded params
savedir='/Volumes/UniWD/AnaNucleus/numPS/';
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
    savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
end
pos=1;
numVars=6;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','k','h','g','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'log \mum^2/ms','log \mum^2/ms','\muM','ms^{-1}','%','\mum','','','ms^{-1}','\muM','ms','ms','ms'};

%% Closest to a sweetspot
tol = 0;
minfdhm=min(dataStruct(cellit).nfdhm)*(1-tol);
maxfdhm=max(dataStruct(cellit).nfdhm)*(1+tol);
minttp=min(dataStruct(cellit).nttpeak)*(1-tol);
maxttp=max(dataStruct(cellit).nttpeak)*(1+tol);
mv=3;
dcfdhm=results(:,12-mv);
dcttp=results(:,11-mv)*10-182;
sweetspot=2*(dcfdhm>minfdhm&dcfdhm<maxfdhm)+3*(dcttp>minttp&dcttp<maxttp);
closeres=find(sweetspot);
possSpots=sum((sweetspot==5))
if mv>0
    result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
    result_units={'log \mum^2/ms','log \mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
end
% resrange=find(results())
for ii=10%[4:8,14:16]
figure
hold on
% plot(dcttp,dcfdhm,results(:,7),'x','MarkerSize',10,'LineWidth',2)
% scatter(results(:,9),results(:,10),10,results(:,7),'LineWidth',2)
scatter(results(:,11-mv)*10-182,results(:,12-mv),10,results(:,ii-mv),'LineWidth',2)
plot([minttp minttp],[minfdhm maxfdhm],'r')
plot([maxttp maxttp],[minfdhm maxfdhm],'r')
plot([minttp maxttp],[maxfdhm maxfdhm],'r')
plot([minttp maxttp],[minfdhm minfdhm],'r')
title(result_labels{ii-mv})
% plot([minttp minttp],[minfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([maxttp maxttp],[minfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([minttp maxttp],[maxfdhm maxfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% plot([minttp maxttp],[minfdhm minfdhm],'Color',[0.8 0 0],10,'LineWidth',2)
% axis([minttp*0.7 maxttp*1.4 minfdhm*0.7 maxfdhm*1.3])
end
%%
figure
plot(results(:,1:5),results(:,12-mv)./results(:,11-mv),'x','LineWidth',2)
%% Effect of changes to params on quals - heatmap
%load('/Users/hhunt1/Documents/huia/parPS_hl1')
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','k','h','g','Max [Ca^{2+}]_{nuc}','TTP','FDHM','Time to base'};
result_units={'log \mum^2/ms','log \mum^2/ms','\muM','ms^{-1}','%','\mum','','','ms^{-1}','\muM','ms','ms','ms'};
save_labels={'D_c','D_b','b_t','K','NPC','a','k','h','g','Max','ttpeak','FDHM','ttbase','hcyt','lcyt'};
numpervar=6;%11;%round(exp(log(size(results,1))/6));
un=cell(1,6);
for unvar=[1:9]
    un{unvar}=unique(results(:,unvar));
end
resrange=1:size(results,1);
for fixD=[1]
    resrangepreloop=intersect(find(results(:,fixD)==un{fixD}(1)),resrange);
end
for qual=11:12
    for ii=7:9
        for jj=setdiff(7:9,ii)
            for fixvar=setdiff(7:9,[ii,jj])
                resrangeinloop=intersect(find(results(:,fixvar)==un{fixvar}(5)),resrangepreloop);
            end
    figure('position',[0,200,560,450])
    if 1
    dc=reshape(results(resrangeinloop,qual),[numpervar,numpervar]);
    x=unique(results(resrangeinloop,ii),'stable');
    y=unique(results(resrangeinloop,jj),'stable');
    surf(x,y,dc','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
    view(90,90)
    caxis([min(results(:,qual)) 850])%max(results(:,qual))])
    axis([min(x) max(x) min(y) max(y)])
%     set(gca,'XScale','log')
    else
    x = results(resrangeinloop,ii);
    y = results(resrangeinloop,jj);
    dc= results(resrangeinloop,qual);
    scatter(x,y,4500,dc,'filled')
    set(gca,'FontSize',baseSize*0.8)
    axis('square')
%     axis([min(x-2) max(x+2) min(y-0.25) max(y+0.25)])
    end
    if jj==9
        set(gca,'YScale','log')
    end
    if ii==9
        set(gca,'XScale','log')
    end
    title(strcat(result_labels{qual}),'FontSize',25)
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',20)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',20)
    colorbar
%     saveas(gcf,strcat(savedir,'chapterfig_cyttrans-',num2str(pos),'_',save_labels{qual},'_','nobuf','.eps'),'epsc')
        end
    end
end
