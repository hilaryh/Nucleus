function [sol,xfdhm,actual_fdhm,actual_ttpeak,actual_nfdhm,actual_nttpeak]=fitDinNuc_spline_num(dataStruct,cellit,a,redoBool)
if (1-redoBool) && ~isempty(dataStruct(cellit).Dvals)
    sol=dataStruct(cellit).Dvals;
    return
end
addpath('/Users/hhunt1/Documents/Nucleus/cutNuc1DData/')
maxcvals=dataStruct(cellit).maxtime;
fdhmcvals=dataStruct(cellit).fdhmtime;
cutangle=dataStruct(cellit).cutangle;
% nucr=(dataStruct(cellit).nucleus(1)-10):(dataStruct(cellit).nucleus(end)+10);
nucr=(dataStruct(cellit).nucleus(1)+7):(dataStruct(cellit).nucleus(end)-7);
% Compare sim and data at each peak. Find D that fits.
ncPos=1;
mccpos=dataStruct(cellit).mccpos{ncPos};
fits=dataStruct(cellit).fstfit{ncPos};
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
actual_nfdhm=x;
actual_nttpeak=x;
actual_fdhm=x;
actual_ttpeak=x;
    for peakit=1:numPeaks
        peakstart = max(round(mccpos(peakit)-0.1*meanpeakDist)+dataStruct(cellit).init,1);
        peakend=min(round(mccpos(peakit)+0.9*meanpeakDist)+dataStruct(cellit).init,size(dataStruct(cellit).smoothedImg,2));
        nucSliceO=dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(nucr),peakstart:peakend);
        
        % Do exactly the same calculation as was done to cytosol...
        cytSlice=median(dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(dataStruct(cellit).nucleus(1:10)-15),peakstart:peakend));
        peakMin=median(cytSlice((end-1000):end));
        peakMax=max(cytSlice);
        nucSlice=(0.9*nucSliceO+0.1*peakMax-peakMin)/(peakMax-peakMin);
        [~,tomaxc]=max(cytSlice);
        tofdhm=range(find(cytSlice>=((max(cytSlice)+min(cytSlice))/2)));
        actual_fdhm(peakit)=tofdhm*0.26;
        actual_ttpeak(peakit)=tomaxc*0.26-fits{peakit,1}.breaks(1);
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
        thisnuc=@(D)callcylfun(D,a,crtpts(peakit,:),fits(peakit,:),...
            tinterval*timeStep,tend*timeStep,nucSlicex);
        %% Fit to curve of D vs time to max
        nucTrans=median(nucSlice(round((end/2-3):(end/2+3)),:));
        [~,tomax]=max(nucTrans);
        x(peakit)=((tomax(1)*timeStep-maxcvals(3))/maxcvals(1))^(1/maxcvals(2));
        if abs(imag(x(peakit)))>0
            warning('Imaginary MAX')
        end
        tofdhm=range(find(nucTrans>=((max(nucTrans)+min(nucTrans))/2)));
        actual_nfdhm(peakit)=tofdhm*0.26;
        actual_nttpeak(peakit)=tomax*0.26-fits{peakit,1}.breaks(1);
        xfdhm(peakit)=((tofdhm(1)*timeStep-fdhmcvals(3))/fdhmcvals(1))^(1/fdhmcvals(2));
        if abs(imag(xfdhm(peakit)))>0
            warning('Imaginary FDHM')
            xfdhm(peakit)=x(peakit);
        end
        strcat('Stageit:',num2str(cellit),'peak:',num2str(peakit),'/',num2str(numPeaks),'D: ',num2str(x(peakit)))
        if 0
            simslice=thisnuc(xfdhm(peakit));
    %         maxq=max(median(matsimslice(:,:)));
    %         rescalefactor=max(max(nucSlice(:,peakstart:10:(peakstart+4000))))/maxq;
            nlength=size(nucSlice,1);
            slength=size(simslice,2);
            svals=[ones(1,10),round(1:(slength/(nlength-20)):slength),slength*ones(1,10)];
            simslice2=num2cell(simslice(:,svals),2)';
            % combine nucSlice and simslice
            ncslice=num2cell(nucSlice(:,1:tinterval:(tend+1))',2)';
            size(ncslice{1});
            size(simslice2{1});
            comslice=cellfun(@(x,y)[x;y],ncslice,simslice2,'Un',0);
    %         x=gamultiobj(leastsquarescomp,1,[],[],[],[],[0],[1e4]);
            makemovie1D(dataStruct(cellit).timeStep,comslice,nucSlicex,cellit,peakit,...
                strcat(dataStruct(cellit).name(1:8),'-com'),10*spaceStep,10*spaceStep,cutangle)
        end
    end
    sol=x;
end
function sol=callcylfun(D,a,crtpts,fits,tint,tmax,nucSlicex)
    if D>0
        [cutca,~,~]=nonlinearDiff(D,[],fits,crtpts,tint,tmax);
        sol=cutca(:,:,1);
        sol=sol(:,[(end:-1:1) (2:end)]);
    else
        peaklength=tmax/tint+1;
        sol=5e4*ones(round(peaklength),size(nucSlicex,2));
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