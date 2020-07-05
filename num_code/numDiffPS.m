%% Parameter sweep of numerical calcium diffusion in cell nucleus
% function numDiffPS()
numEachVar=2;
numSteps = numEachVar^6;
ns = numEachVar-1;
param_type=1;
if param_type
    % Choose parameters on grid
    D_u = 10.^(6*repmat(0:(1/ns):1,[1,numEachVar^6])-4);
    D_b = 10.^(6*repmat(reshape(repmat((0:(1/ns):1),[numEachVar^4,1]),numEachVar^5,1),[numEachVar,1])-4);
    b_t = 1e4*reshape(repmat(0:(1/ns):1,[numEachVar^5,1]),numEachVar^6,1);
    K   = 1e4*repmat(reshape(repmat((0:(1/ns):1),[numEachVar^3,1]),numEachVar^4,1),[numEachVar^2,1]);
    NPC = 2*repmat(reshape(repmat((0:(1/ns):1),[numEachVar^2,1]),numEachVar^3,1),[numEachVar^3,1]);
    a   = 2+10*repmat(reshape(repmat((0:(1/ns):1),[numEachVar,1]),numEachVar^2,1),[numEachVar^4,1]);

else
% Choose quasi-random parameters 
rng default
p = sobolset(6,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
p = net(p,numSteps);

D_u = 10.^(6*p(:,1)-4);
D_b = 10.^(6*p(:,2)-4);
b_t = 1000*p(:,3);
K   = 1e4*p(:,4);
NPC = 2*p(:,5);
a   = 10*p(:,6)+2;
end
savedir='/Users/hhunt1/Documents/Nucleus/cutNuc1DData/num_plots/';
if ~exist(savedir,'dir')
    savedir='/home/hhunt1/AnaNuc/1Dsims/numerical/plots/';
end

params=[D_b,b_t,K,NPC,a];
for pos = 1%1:4
fits=dataStruct(cellit).fstfit{3};
crtpts=dataStruct(cellit).crtpts{3};
tint=10;
tmax=3000;
results=zeros(numSteps,10);
for ii=1:numSteps
    [sol,r,t]       =nonlinearDiff(D_u(ii),params(ii,:),fits,crtpts,tint,tmax);
    
    nucDbl          = sol(:,1,1);
    [maxND,ttpeak]  = max(nucDbl);
    fdhm            = range(find(nucDbl>((maxND+min(nucDbl(:)))/2)))*tint;
    fdnb            = range(find(nucDbl>(0.1*maxND+0.9*min(nucDbl(:)))))*tint;
    if isempty(fdhm)
        fdhm=0;
    end
    if isempty(fdnb)
        fdnb=0;
    end
try
    results(ii,:)   = [D_u(ii),params(ii,:),maxND,ttpeak,fdhm,fdnb];
catch
    warning('Something wrong when recording results')
end
end
if 0
%% Create figures of results
label='nc1';
numVars=6;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','Time to peak','FDHM','Time to base'};
result_units={'log \mum^2/ms','log \mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
for ii=1:(numVars-1)
    for jj=(ii+1):numVars
for qual = numVars+(1:4)
    x = results(:,ii);
    y = results(:,jj);
    dc= results(:,qual);
    if ii==1 || ii==2
        x = log(x);
    end
    if jj==2
        y = log(y);
    end
    f = fit([x,y],dc,'linearinterp');
    
    figure('position',[0,200,560,450])
    plot(f,'Style','contour')
    set(gca,'FontSize',16)
    view(90,-90)
    axis('square')
    colorbar
    title(strcat(result_labels{qual}),'FontSize',25)
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',20)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',20)
      saveas(gcf,strcat(savedir,label,'contour_pos-',num2str(pos),'_',num2str(ii),'_',num2str(jj),'_',...
    num2str(qual),'.eps'),'epsc')
saveas(gcf,strcat(savedir,label,'contour_pos-',num2str(pos),'_',num2str(ii),'_',num2str(jj),'_',...
    num2str(qual),'.fig'),'epsc')
end
    end
end
end
save(strcat(savedir,'PS',num2str(pos)),'results')
end

%% Create figures of results
numVars=6;
result_labels={'D_c','D_b','b_{tot}','K_b','NE','a','Max [Ca^{2+}]_{nuc}','Time to peak','FDHM','Time to base'};
result_units={'log \mum^2/ms','log \mum^2/ms','\muM','ms^{-1}','%','\mum','\muM','ms','ms','ms'};
for ii=1:(numVars-1)
    for jj=(ii+1):numVars
for qual = numVars+(1:4)
    x = results(:,ii);
    y = results(:,jj);
    dc= results(:,qual);
    if ii==1 || ii==2
        x = log(x);
    end
    if jj==2
        y = log(y);
    end
    
    figure('position',[0,200,560,450])
    scatter3(x,y,dc,50,dc,'filled')
    set(gca,'FontSize',16)
    view(90,-90)
    axis('square')
    colorbar
    title(strcat(result_labels{qual}),'FontSize',25)
    ylabel(strcat(result_labels(jj),{' '},'(',result_units(jj),')'),'FontSize',20)
    xlabel(strcat(result_labels(ii),{' '},'(',result_units(ii),')'),'FontSize',20)
      saveas(gcf,strcat(savedir,'scatter_pos-',num2str(pos),'_',num2str(ii),'_',num2str(jj),'_',...
    num2str(qual),'.eps'),'epsc')
saveas(gcf,strcat(savedir,'scatter_pos-',num2str(pos),'_',num2str(ii),'_',num2str(jj),'_',...
    num2str(qual),'.fig'),'epsc')
end
    end
end

