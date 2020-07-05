function sol=fitDinNuc_spline2(dataStruct,cellit,a,l,redoBool)
if (1-redoBool) && ~isempty(dataStruct(cellit).Dvals)
    sol=dataStruct(cellit).Dvals;
    return
end
maxcvals=dataStruct(cellit).maxtime;
fdhmcvals=dataStruct(cellit).fdhmtime;
cutangle=dataStruct(cellit).cutangle;
% nucr=(dataStruct(cellit).nucleus(1)-10):(dataStruct(cellit).nucleus(end)+10);
nucr=(dataStruct(cellit).nucleus(1)+7):(dataStruct(cellit).nucleus(end)-7);
% Compare sim and data at each peak. Find D that fits.
ncPos=3;
mccpos=dataStruct(cellit).mccpos{ncPos};
fits=dataStruct(cellit).fits{ncPos};
crtpts=dataStruct(cellit).crtpts{ncPos};
numPeaks=size(mccpos,2);
meanpeakDist=mean(diff(mccpos));
if isempty(meanpeakDist)
    warning('No peaks? Check fitD line 10')
    sol=0;
    return
end
timeStep=dataStruct(1).timeStep;
spaceStep=dataStruct(1).spaceStep;
x=zeros(1,numPeaks);
xfdhm=x;
    for peakit=1:numPeaks
        peakstart = max(round(mccpos(peakit)-0.1*meanpeakDist)+dataStruct(cellit).init,1);
        peakend=min(round(mccpos(peakit)+0.9*meanpeakDist)+dataStruct(cellit).init,size(dataStruct(cellit).smoothedImg,2));
        nucSliceO=dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(nucr),peakstart:peakend);
        
        % Do exactly the same calculation as was done to cytosol...
        cytSlice=median(dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(dataStruct(cellit).nucleus(1:10)-15),peakstart:peakend));
        peakMin=median(cytSlice((end-1000):end));
        peakMax=max(cytSlice);
        nucSlice=(0.9*nucSliceO+0.1*peakMax-peakMin)/(peakMax-peakMin);
        
        % Ljubojevic F/F0 to uM equation
%         minnucF=min(median(nucSliceO(21:(end-20),:)',2));
%         R_f=max(median(nucSliceO(21:(end-20),:)',2))/minnucF;
%         K_d=2.200;
%         ca_rest=0.200;
%         nucSlice2=K_d*(nucSliceO/minnucF*(K_d+R_f*ca_rest)-(K_d+ca_rest))./(R_f*(K_d+ca_rest)-nucSliceO/minnucF*(K_d+R_f*ca_rest));
        % simulated nuc should come back with the same time and space
        % intervals
        tinterval=15;
        tend=6000;
        nucSlicex=(ceil(-size(nucr,2)/2):floor(size(nucr,2)/2))*spaceStep;
%         thisnuc=@(D)callcylfun(D,l,a,cutangle,crtpts(peakit,:),fits(peakit,:),0,...
%             tinterval*timeStep,(peakend-peakstart)*timeStep,nucSlicex);
        thisnuc=@(D)callcylfun(D,l,a,cutangle,crtpts(peakit,:),fits(peakit,:),0,...
            tinterval*timeStep,tend*timeStep,nucSlicex);
        if 0
        %% Matlab optimiser method
%         fminsearch(
% Need to line up ts and then least squares each timestep that matches.
%         leastsquarescomp=@(D)sum(cell2mat(cellfun(@(x,y)(x(lcut:hcut)-y(lcut:hcut)).^2,num2cell(nucSlice(:,peakstart:peakend)',1),...
%             cellfun(@(x)x*rescalefactor,thisnuc(D),'Un',0),'Un',0)));
%         x(peakit)=fminsearch(leastsquarescomp,12873); %33053 %19261 %15544
        leastsquarescompfct=@(D)leastsquarescomp(D,rescalefactor,thisnuc,...
            nucSlice,peakstart,peakend,lcut,hcut);
        x(peakit)=fminsearch(leastsquarescompfct,12873);
        else
        %% Fit to curve of D vs time to max
        nucTrans=median(nucSlice(round((end/2-3):(end/2+3)),:));
        [~,tomax]=max(nucTrans);
        x(peakit)=((tomax(1)*timeStep-maxcvals(3))/maxcvals(1))^(1/maxcvals(2));
        tofdhm=range(find(nucTrans>=((max(nucTrans)+min(nucTrans))/2)));
        xfdhm(peakit)=((tofdhm(1)*timeStep-fdhmcvals(3))/fdhmcvals(1))^(1/fdhmcvals(2));
        end
        strcat('Stageit:',num2str(cellit),'peak:',num2str(peakit),'/',num2str(numPeaks),'D: ',num2str(x(peakit)))
        simslice=thisnuc(xfdhm(peakit));
        matsimslice=cell2mat(simslice')';
%         maxq=max(median(matsimslice(:,:)));
%         rescalefactor=max(max(nucSlice(:,peakstart:10:(peakstart+4000))))/maxq;
        nlength=size(nucSlice,1);
        slength=size(simslice{1},2);
        svals=[ones(1,10),round(1:(slength/(nlength-20)):slength),slength*ones(1,10)];
        simslice2=cellfun(@(x)x(svals),simslice,'Un',0);
        % combine nucSlice and simslice
        ncslice=num2cell(nucSlice(:,1:10:4001)',2)';
        size(ncslice{1})
        size(simslice2{1})
        comslice=cellfun(@(x,y)[x;y],ncslice,simslice2,'Un',0);
%         x=gamultiobj(leastsquarescomp,1,[],[],[],[],[0],[1e4]);
        makemovie1D(dataStruct(cellit).timeStep,comslice,nucSlicex,cellit,peakit,...
            strcat(dataStruct(cellit).name(1:8),'-com'),10*spaceStep,10*spaceStep,cutangle)
    end
    sol=x;
end
function sol=callcylfun(D,l,a,cutangle,crtpts,fits,varchoice,tint,tmax,nucSlicex)
    if D>0
        [cutca,h,~,~,~,~]=sim_nucleus_cyl_wz_spline(D,l,a,cutangle,crtpts,fits,tint,tmax,1);
        if varchoice==0
            sol=cellfun(@(x)x,cutca,'UniformOutput',false);
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