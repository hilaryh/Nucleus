function sol=fitD(dataStruct,stageit,cutangle,nucSlice,nucSlicex,hlen,llen,a,l,typeit,cellit,redoBool)
if (1-redoBool) && ~isempty(dataStruct(stageit).Dvals)
    sol=dataStruct(stageit).Dvals;
    return
end
% Compare sim and data at each peak. Find D that fits.
numPeaks=size(dataStruct(stageit).mccpos,2);
mccpos=dataStruct(stageit).mccpos;
meanpeakDist=mean(diff(mccpos));
if isempty(meanpeakDist)
    warning('No peaks? Check fitD line 10')
    sol=0;
    return
end
caconclen=size(nucSlice,2);
timeStep=dataStruct(1).timeStep;
spaceStep=dataStruct(1).spaceStep;
qvars=dataStruct(stageit).qvars;
lcut=find(nucSlicex<=llen,1,'last')+5;
hcut=find(nucSlicex>=(max(nucSlicex)-hlen),1,'first')-5;
x=zeros(1,numPeaks);
    for peakit=1:numPeaks
        maxq=max(q((0:meanpeakDist)*timeStep/1000,qvars(peakit,2),qvars(peakit,1),...
            qvars(peakit,3),qvars(peakit,4),0));
        peakstart = max(round(mccpos(peakit)-0.1*meanpeakDist),1);
        peakend=min(round(mccpos(peakit)+0.9*meanpeakDist),size(nucSlice,2));
        rescalefactor=max(max(nucSlice(:,peakstart:peakend)))/maxq;
        % simulated nuc should come back with the same time and space
        % intervals
        thisnuc=@(D)callcylfun(D,l,a,cutangle,dataStruct(stageit).qvars(peakit,:),0,...
            timeStep/1000,(peakend-peakstart)*timeStep/1000,nucSlicex);
%         fminsearch(
% Need to line up ts and then least squares each timestep that matches.
%         leastsquarescomp=@(D)sum(cell2mat(cellfun(@(x,y)(x(lcut:hcut)-y(lcut:hcut)).^2,num2cell(nucSlice(:,peakstart:peakend)',1),...
%             cellfun(@(x)x*rescalefactor,thisnuc(D),'Un',0),'Un',0)));
%         x(peakit)=fminsearch(leastsquarescomp,12873); %33053 %19261 %15544
        leastsquarescompfct=@(D)leastsquarescomp(D,rescalefactor,thisnuc,...
            nucSlice,peakstart,peakend,lcut,hcut);
        x(peakit)=fminsearch(leastsquarescompfct,12873);
        strcat('Stageit:',num2str(stageit),'peak:',num2str(peakit),'/',num2str(numPeaks),'D: ',num2str(x(peakit)))
        simslice=thisnuc(x(peakit));
        simslice2=cellfun(@(x)x*rescalefactor,simslice,'Un',0);
        % combine nucSlice and simslice
        comslice=cellfun(@(x,y)[x;y],num2cell(nucSlice(:,peakstart:peakend)',2)',simslice2,'Un',0);
%         x=gamultiobj(leastsquarescomp,1,[],[],[],[],[0],[1e4]);
        makemovie1D(dataStruct(stageit).timeStep,comslice,nucSlicex,cellit,peakit,...
            strcat(dataStruct(stageit).name(1:8),'-com'),hlen,llen,cutangle)
    end
    sol=x;
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

function sol= leastsquarescomp(D,rescalefactor,thisnuc,...
            nucSlice,peakstart,peakend,lcut,hcut)
% cf1 and ca should be cell arrays of length peakend-peakstart+1 with each
% cell being a slice of the nucleus at each time point.
cf1=cellfun(@(x)x*rescalefactor,thisnuc(D),'Un',0);
ca=num2cell(nucSlice(:,peakstart:peakend)',2)';
if sum(size(cf1)==size(ca))~=2
    warning(strcat('Something wrong with sizes - cf1: ',num2str(size(cf1,1)),...
        ',',num2str(size(cf1,2)),' ca: ',num2str(size(ca,1)),...
        ',',num2str(size(ca,2))))
end
cf2=cellfun(@(x,y)(x(lcut:hcut)-y(lcut:hcut)).^2,ca,cf1,'Un',0);
sol=sum(cell2mat(cf2));
end