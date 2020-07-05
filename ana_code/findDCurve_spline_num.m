function [maxsol,fdhmsol,Drange,maxtime]=findDCurve_spline_num(dataStruct,cellit,a,redoBool)
    cutangle=dataStruct(cellit).cutangle;
    timeStep=dataStruct(cellit).timeStep;
    nucSlicex=dataStruct(cellit).nucSlicex;
    tint=50*timeStep;
    ncPos=1;
    thisnuc=@(D)callcylfun(D,a,mode(dataStruct(cellit).crtpts{ncPos},1),dataStruct(cellit).fstfit{ncPos},...
            tint,11000*timeStep,nucSlicex);
    numpts=10;
    Drange=logspace(-2.5,2,numpts);
    maxtime=zeros(size(Drange));
    fdhm=maxtime;
    nucii=cell(size(Drange));
    for ii=1:numpts
        nucii{ii}=thisnuc(Drange(ii));
        nucDbl=nucii{ii}(:,1);
        [maxND,maxtime(ii)]=max(nucDbl);
        fdhm(ii)=range(find(nucDbl>((maxND+min(nucDbl(:)))/2)))*tint;
        if maxtime(ii)<0
            warning('Something wrong with finding max')
        end
    end
    tomaxfit=fit(Drange',tint*maxtime','power2');
    maxsol=coeffvalues(tomaxfit);
    tofdhmfit=fit(Drange',fdhm','power2');
    fdhmsol=coeffvalues(tofdhmfit);
    t1=logspace(-3,3,1e3);
    t2=tomaxfit.a*t1.^tomaxfit.b+tomaxfit.c;
    figure
    semilogx(t1,t2,'LineWidth',2)
    hold on
    semilogx(Drange,tint*maxtime,'x','MarkerSize',20,'LineWidth',2)
    ylabel('Time to peak in nucleus (ms)','FontSize',16)
    xlabel('Diffusion coefficient (\mum^2ms^{-1})','FontSize',16)
    legend({'fit','simulation'},'FontSize',16)
    if exist(dataStruct(cellit).fpath(1:(end-1)),'dir')
    saveas(gcf,strcat(dataStruct(cellit).fpath,'graphs/findDCurve-',...
        dataStruct(cellit).name,'.png'))
    end
    t1=logspace(-3,3,1e3);
    t2=tofdhmfit.a*t1.^tofdhmfit.b+tofdhmfit.c;
    figure
    semilogx(t1,t2,'LineWidth',2)
    hold on
    semilogx(Drange,fdhm,'x','MarkerSize',20,'LineWidth',2)
    ylabel('FDHM in nucleus (ms)','FontSize',16)
    xlabel('Diffusion coefficient (\mum^2ms^{-1})','FontSize',16)
    legend({'fit','simulation'},'FontSize',16)
    if exist(dataStruct(cellit).fpath(1:(end-1)),'dir')
    saveas(gcf,strcat(dataStruct(cellit).fpath,'graphs/findDCurve_fdhm-',...
        dataStruct(cellit).name,'.png'))
    end
    
    % -> sol(D) = cval(1)*D^cval(2)+cval(3)
end
function sol=callcylfun(D,a,crtpts,fits,tint,tmax,nucSlicex)
    if D>0
        [cutca,~,~]=nonlinearDiff(D,[],fits,crtpts,tint,tmax);
        sol=cutca(:,:,1);
    else
        sol=5e4*ones(round(tmax/tint)+1,size(nucSlicex,2));
    end
end