function [maxsol,fdhmsol,Drange,maxtime]=findDCurve_spline(dataStruct,stageit,a,l,typeit,redoBool)
    cutangle=dataStruct(stageit).cutangle;
    timeStep=dataStruct(stageit).timeStep;
    nucSlicex=dataStruct(stageit).nucSlicex;
    tint=50*timeStep;
    ncPos=3;
    thisnuc=@(D)callcylfun(D,l,a,cutangle,mode(dataStruct(stageit).crtpts{ncPos}),dataStruct(stageit).fits{ncPos},0,...
            tint,11000*timeStep,nucSlicex);
    numpts=5;
    Drange=logspace(-2.5,0,numpts);
    maxtime=zeros(size(Drange));
    fdhm=maxtime;
    nucii=cell(size(Drange));
    for ii=1:numpts
        nucii{ii}=thisnuc(Drange(ii));
        nucDbl=cell2mat(nucii{ii}');
        [maxND,maxtime(ii)]=max((nucDbl));
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
    if exist(dataStruct(stageit).fpath(1:(end-1)),'dir')
    saveas(gcf,strcat(dataStruct(stageit).fpath,'graphs/findDCurve-',...
        dataStruct(stageit).name,'.png'))
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
    if exist(dataStruct(stageit).fpath(1:(end-1)),'dir')
    saveas(gcf,strcat(dataStruct(stageit).fpath,'graphs/findDCurve_fdhm-',...
        dataStruct(stageit).name,'.png'))
    end
    
    % -> sol(D) = cval(1)*D^cval(2)+cval(3)
end
function sol=callcylfun(D,l,a,cutangle,crtpts,fits,varchoice,tint,tmax,nucSlicex)
    if D>0
        [cutca,h,~,~,~,~]=sim_nucleus_cyl_wz_spline(D,l,a,cutangle,crtpts,fits,tint,tmax,1);
        if varchoice==0
%             cutcax=[-h(end:-1:2) h];
%             newNSx=round(nucSlicex-max(nucSlicex)/2,2);
%             [~,idx]=ismember(newNSx,round(cutcax,2));
%             idx=nonzeros(idx);
%             sol=cellfun(@(x)x(idx),cutca,'Un',0);
            sol=cellfun(@(x)median(x(800:1009)),cutca,'UniformOutput',false);
%             sol=[centreAv{:}];
        else
            sol=h;
        end
    else
%         warning(strcat('D=',num2str(D)))
        peaklength=tmax/tint+1;
        sol=cell(1,round(peaklength));
        for ii=1:size(sol,2)
            sol{ii}=5e4*ones(1,size(nucSlicex,2));
        end
    end
end