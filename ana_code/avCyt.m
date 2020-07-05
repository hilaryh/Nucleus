% H Hunt 2019
% Called by cutNuc1DData.m
% Mba=qvars=[beat,M,alpha,c_0,padding_length]
function [Mba,mccpos,plusff,thisPeak]=avCyt(dataStruct,stageit)
    timeSteps=size(dataStruct(stageit).img,2)-dataStruct(stageit).init-1;
    plusff=[];
    cytosol=dataStruct(stageit).cytosol;
%     cytBd=boundary(mod(cytosol,sqSize)+1,ceil(cytosol/sqSize),1);
    % Find peaks
    cytca=median(dataStruct(stageit).smoothedImg(dataStruct(stageit).xrng(cytosol),(dataStruct(stageit).init+1):end));
    [fpks,fpos,~,fprom]=findpeaks(-cytca,'MinPeakDistance',2800);
    [~,promorder]=sort(fprom);
    minProm=fprom(promorder(end-4));
    [~,mccpos]=findpeaks(cytca,'MinPeakDistance',2800,'MinPeakProminence',minProm);
    meanpeakDist=mean(diff(mccpos));
%     % Convert to F/F0?
    c_0=median(cytca([max(round(mccpos-0.1*meanpeakDist),1), max(round(mccpos-0.2*meanpeakDist),1)]));
%     cytca=cytca/c_0;
%     [cytcaco,plusff]=FF0touM(cytca,1,0);
    cytcaco=cytca;
%     c_0=1e3*median(cytcaco([max(round(mccpos-0.1*meanpeakDist),1), max(round(mccpos-0.2*meanpeakDist),1)]));
    numPeaks=size(mccpos,2);
    thisPeak=cell(1,numPeaks);
    Mba=zeros(numPeaks,4);
    Mba(:,4)=c_0*ones(numPeaks,1);
    for peakit=1:numPeaks
        thisPeak{peakit}=median(dataStruct(stageit).smoothedImg(dataStruct(stageit).xrng(cytosol),...
                (max(round(mccpos(peakit)-0.1*meanpeakDist)+dataStruct(stageit).init,1):...
            min(round(mccpos(peakit)+0.9*meanpeakDist)+dataStruct(stageit).init,size(dataStruct(stageit).smoothedImg,2)))));
        [~,diffmax]=max(diff(thisPeak{peakit}));
        flatPart=find(diff(thisPeak{peakit}(1:diffmax))<mean(diff(thisPeak{peakit})),1,'last')-110;
        thisPeak{peakit}=thisPeak{peakit}(flatPart:end);
        tPointToS=dataStruct(stageit).timeStep/1000;
        % Compare widths and amplitudes
        dfdhm=findfdhm(thisPeak{peakit});
        dfd8m=findfd8m(thisPeak{peakit});
        [~,mdatapos]=max(thisPeak{peakit});
        timevec=(0:(size(thisPeak{peakit},2)-1))*tPointToS;
        timevecopt=timevec;%(0:1000)*tPointToS;
        fdhmcomp=@(optvar)abs(dfdhm-findfdhm(q(timevecopt,optvar(1),optvar(2),optvar(3),c_0,0)));
        fd8mcomp=@(optvar)abs(dfd8m-findfd8m(q(timevecopt,optvar(1),optvar(2),optvar(3),c_0,0)));
        basecomp=@(optvar)abs(min(thisPeak{peakit})-min(q(timevecopt,optvar(1),optvar(2),optvar(3),c_0,0)));
        maxcomp=@(optvar)abs(max(thisPeak{peakit})-max(q(timevecopt,optvar(1),optvar(2),optvar(3),c_0,0)));
%         ttpeakcomp=@(optvar)abs(mdatapos-find(thisq==max(thisq),1,'first'));
        newfunqopt=@(optvar)fdhmcomp(optvar)+fd8mcomp(optvar)+maxcomp(optvar)+basecomp(optvar)+1e-2*ttpeakcomp(optvar,c_0,timevecopt,mdatapos);
        optqparams=fminsearch(newfunqopt,[2.4,0.4,0.5]);
        optqparams=fminsearch(newfunqopt,optqparams);
        Mba(peakit,1:3)=optqparams;
        % Align peaks
        thisq=q(timevec,Mba(peakit,1),Mba(peakit,2),Mba(peakit,3),c_0,0);
        [~,mqpos]=max(thisq);
        Mba(peakit,5)=max(mdatapos-mqpos,0);
        if 1
            figure
            plot(timevec*1e3,thisPeak{peakit})
            hold on
            plot(timevec*1e3,thisq)
        end
        
    end
end
function sol=ttpeakcomp(optvar,c_0,timevecopt,mdatapos)
    thisq=q(timevecopt,optvar(1),optvar(2),optvar(3),c_0,0);
    sol=abs(mdatapos-find(thisq==max(thisq),1,'first'));
end
function sol=ttpeakcomp2(optvar,c_0,timevecopt,mdatapos)
    thisq=q2(timevecopt,optvar(1),optvar(2),optvar(3),c_0,0);
    sol=abs(mdatapos-find(thisq==max(thisq),1,'first'));
end

function sol=findfdhm(catrans)
maxcat=max(catrans);
mincat=min(catrans);
% interpolate then find
longcatrans=interp1(1:size(catrans,2),catrans,1:1e-3:size(catrans,2));
fdhm=find(longcatrans>=((maxcat+mincat)/2));
sol=1e-3*(fdhm(end)-fdhm(1));
end
function sol=findfd8m(catrans)
maxcat=max(catrans);
mincat=min(catrans);
% interpolate then find
longcatrans=interp1(1:size(catrans,2),catrans,1:1e-3:size(catrans,2));
fdhm=find(longcatrans>=((maxcat-mincat)*0.8+mincat));
sol=1e-3*(fdhm(end)-fdhm(1));
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