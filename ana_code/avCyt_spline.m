% H Hunt 2019
% Called by cutNuc1DData.m
% Mba=qvars=[beat,M,alpha,c_0,padding_length]
function [crtpts,mccpos,fstfit,thisPeak]=avCyt_spline(dataStruct,cellit,redoBool)
if redoBool==0 && ~isempty(dataStruct(cellit).crtpts)
    crtpts={dataStruct(cellit).crtpts,dataStruct(cellit).crtptsr};
    fstfit={dataStruct(cellit).fits,dataStruct(cellit).fitsr};
    thisPeak={dataStruct(cellit).thisPeak,dataStruct(cellit).thisPeakr};
    return
end
    crtpts=cell(1,2);
    mccpos=crtpts;
    fstfit=crtpts;
    thisPeak=crtpts;
    timeStep=dataStruct(cellit).timeStep;
    % Look at cytosol near nucleus
    cytosol=cell(1,4);
    % Just outside the nuc l
    cytosol{1}=dataStruct(cellit).nucleus(1:10)-15;
    % Just outside the nuc r
    cytosol{2}=dataStruct(cellit).nucleus((end-10):end)+15;
    % Look at space where nucleus transient transforms (Is this where NE
    % actually is?)
    spaceNuc=dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(dataStruct(cellit).nucleus),dataStruct(cellit).init+7240);
    [~,maxl]=max(spaceNuc(1:15));
    [~,maxr]=max(spaceNuc((end-14):end));
    % Just inside the nuc l
    cytosol{3}=dataStruct(cellit).nucleus(maxl);
    % Just inside the nuc r
    cytosol{4}=dataStruct(cellit).nucleus(end-13+maxr);
    % Centre of nucleus
    cytosol{5}=dataStruct(cellit).nucleus(round(end/2));
    for lr=1:size(cytosol,2)
%     cytBd=boundary(mod(cytosol,sqSize)+1,ceil(cytosol/sqSize),1);
    % Find peaks
    cytca=median(dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(cytosol{lr}),(dataStruct(cellit).init+1):end),1);
    try
        [~,~,~,fprom]=findpeaks(-cytca,'MinPeakDistance',2800);
    catch
        
        warning('empty cytca')
    end
    [~,promorder]=sort(fprom);
    minProm=fprom(promorder(end-4));
    [~,mccposit]=findpeaks(cytca,'MinPeakDistance',2800,'MinPeakProminence',minProm);
    meanpeakDist=mean(diff(mccposit));
%     % Convert to F/F0?

    numPeaks=size(mccposit,2);
    thisPeakit=cell(1,numPeaks);
    crtptsit=zeros(numPeaks,14);
    fstfitit=cell(numPeaks,2);
    for peakit=1:numPeaks
        thisPeakit{peakit}=median(dataStruct(cellit).smoothedImg(dataStruct(cellit).xrng(cytosol{lr}),...
                (max(round(mccposit(peakit)-0.1*meanpeakDist)+dataStruct(cellit).init,1):...
            min(round(mccposit(peakit)+0.9*meanpeakDist)+dataStruct(cellit).init,size(dataStruct(cellit).smoothedImg,2)))),1);
        peakMin=thisPeakit{peakit}(end);
        peakMax=max(thisPeakit{peakit});
        thisPeakit{peakit}=(0.9*thisPeakit{peakit}+0.1*peakMax-peakMin)/(peakMax-peakMin);
        % Find critical points/breaks for spline fits
        if max(diff(thisPeakit{peakit}))>1e-3
            % Approx. start of rise
            cpts(1)=find(diff(thisPeakit{peakit})>1e-3,1,'first');
            % Mid rise
            [~,cpts(2)]=max(diff(thisPeakit{peakit}(1:3000)));
            % Peak
            [~,cpts(3)] = max(thisPeakit{peakit}(1:3000));
            % Mid fall
            [~,cpts(4)]=min(diff(thisPeakit{peakit}));
            % end
            cpts(5)=size(thisPeakit{peakit},2);
            
            try
                % Right before rise
                crtptsit(peakit,1)=find(diff(thisPeakit{peakit}(1:cpts(2)))<median(abs(diff(thisPeakit{peakit}))),1,'last')-15;
            catch
                crtptsit(peakit,1)=1;
            end
            if crtptsit(peakit,1)<1
                crtptsit(peakit,1)=1;
            end
            % Right after rise
            crtptsit(peakit,2)=crtptsit(peakit,1)+15;
            % Add cpts and pts between cpts
            crtptsit(peakit,3:(end-4))=round(sort([cpts,cpts(2:(end-1))-diff(cpts(1:(end-1)))/2]));
            % Add in last cpts
            crtptsit(peakit,end-3)=cpts(end);
            % pt between last cpts and interp point
            crtptsit(peakit,end-4)=crtptsit(peakit,end-5)+400;
            ii=size(crtptsit,2)-4;
            f2=thisPeakit{peakit}(crtptsit(peakit,ii):crtptsit(peakit,ii+1));
            f1=(crtptsit(peakit,ii):(crtptsit(peakit,ii+1)))*timeStep;
            if isempty(f1)&&peakit==numPeaks
                crtptsit(end,:)=[];
                fstfitit(end,:)=[];
                numPeaks=numPeaks-1;
            else
                
            try 
                expfit=fit(f1',f2','exp2');
            catch
                 expfit=fit(f1',f2','poly3');
            end
            crtptsit(peakit,end-2:end)=[(crtptsit(peakit,end-4)+crtptsit(peakit,end-3))/2, ...
                (crtptsit(peakit,end-4)+3*crtptsit(peakit,end-3))/4, ...
                (3*crtptsit(peakit,end-4)+crtptsit(peakit,end-3))/4];
            crtptsit(peakit,:)=sort(crtptsit(peakit,:));
            xpts=crtptsit(peakit,1:(end-4));
            try
                f = fit(crtptsit(peakit,:)'*timeStep,[thisPeakit{peakit}(xpts) expfit(timeStep*crtptsit(peakit,end-3:end))']','smoothingspline','SmoothingParam',0.07);
            catch
                warning(strcat('peakit:',num2str(peakit),'xpts:',num2str(size(xpts)),'crtptsit:',num2str(size(crtptsit))))
                warning(strcat('xpts:',num2str(xpts)))
                warning(strcat('thisPeak:',num2str(size(thisPeak{peakit}))))
            end
            fstfitit{peakit,1}=coeffvalues(f);
            
            
            
            if 0
                cy = [240,228,66]/255;
                cb = [0 114 178]/255;
                cgold= [230 159 0]/255;
                timevec=(0:(size(thisPeakit{peakit},2)-1))*timeStep;
                figure
                plot(timevec,thisPeakit{peakit},'LineWidth',5,'Color',cb)
                hold on
                plot((crtptsit(peakit,1):crtptsit(peakit,end))*timeStep,f((crtptsit(peakit,1):crtptsit(peakit,end))*timeStep),'LineWidth',2,'Color',cgold)
                plot(crtptsit(peakit,1:(end))*timeStep,thisPeakit{peakit}(round(crtptsit(peakit,1:(end)))),'x','LineWidth',2,'MarkerSize',10,'Color',co)
%                 plot((crtptsit(peakit,(end-1)):crtptsit(peakit,end))*timeStep,fstfitit{peakit,2}(round((crtptsit(peakit,(end-1)):crtptsit(peakit,end))*timeStep)))
                if exist('/home/hhunt1/AnaNuc/1Dsims/data/graphs','dir')
                    saveas(gcf,strcat('/home/hhunt1/AnaNuc/1Dsims/data/graphs/',...
                        'fit-',num2str(cellit),'-',num2str(peakit),'.png'))
                end
            end
            end
        end
    end
for peakit=1:numPeaks
if crtptsit(peakit,1)==0
crtptsit(peakit,:)=[];
fstfitit(peakit,:)=[];
thisPeakit(peakit)=[];
end
end
crtptsit=crtptsit*timeStep;
crtpts{lr}=crtptsit;
mccpos{lr}=mccposit;
fstfit{lr}=fstfitit;
thisPeak{lr}=thisPeakit;
    end
end

function [sol,plusff]=FF0touM(img,cyt,plusff)
    % If img range isn't within certain bounds, the equation doesn't work.
    % Think about this further.
if plusff~=0
    nimg=img+plusff;
        img=nimg/min(nimg);
else
    if max(img/min(img))>5
        imgo=img;
        maxff0=2.8;
        plusff=(max(img)-maxff0*min(img))/(maxff0-1);
        nimg=img+plusff;
        img=nimg/min(nimg);
    else
        img=img/min(img);
    end
end
    if cyt
        K_d_cyt = 1.1e-3;
        R_f_cyt= 7;
       sol = K_d_cyt*(img*(K_d_cyt+R_f_cyt*0.0001)-(K_d_cyt+0.0001))./(R_f_cyt*(K_d_cyt+0.0001)-img*(K_d_cyt+R_f_cyt*0.0001));
    else
        K_d_nuc = 1.26e-3;
        R_f_nuc=13;
       sol = K_d_nuc*(img*(K_d_nuc+R_f_nuc*0.0001)-(K_d_nuc+0.0001))./(R_f_nuc*(K_d_nuc+0.0001)-img*(K_d_nuc+R_f_nuc*0.0001));
    end
end
