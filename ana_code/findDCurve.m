function [sol,Drange,maxtime]=findDCurve(dataStruct,stageit,cutangle,nucSlicex,a,l,typeit,redoBool)
    timeStep=dataStruct(stageit).timeStep;
    thisnuc=@(D)callcylfun(D,l,a,cutangle,mode(dataStruct(stageit).qvars),0,...
            timeStep/1000,11*timeStep,nucSlicex);
    numpts=5;
    Drange=logspace(-2,4,numpts);
    maxtime=zeros(size(Drange));
    nucii=cell(size(Drange));
    for ii=1:numpts
        nucii{ii}=thisnuc(Drange(ii));
        nucDbl=cell2mat(nucii{ii}');
        [~,maxtime(ii)]=max(median(nucDbl(:,(end/2-5):(end/2+5)),2));
    end
    tomaxfit=fit(Drange',maxtime'*timeStep,'power2');
    sol=coeffvalues(tomaxfit);
    figure
    plot(tomaxfit,Drange,maxtime*timeStep)
    xlabel('D')
    ylabel('Time to peak in nucleus (ms)')
    saveas(gcf,strcat(dataStruct(stageit).fpath,'graphs/findDCurve-',...
        dataStruct(stageit).name,'.png'))
    
    % -> sol(D) = cval(1)*D^cval(2)+cval(3)
end
function sol=callcylfun(D,l,a,cutangle,qvars,varchoice,tint,tmax,nucSlicex)
    if D>0
        [cutca,h]=sim_nucleus_cyl(D,l,a,cutangle,qvars,tint,tmax);
        if varchoice==0
            cutcax=[-h(end:-1:2) h];
            newNSx=round(nucSlicex-max(nucSlicex)/2,2);
            [~,idx]=ismember(newNSx,round(cutcax,2));
            idx=nonzeros(idx);
            sol=cellfun(@(x)x(idx),cutca,'Un',0);
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