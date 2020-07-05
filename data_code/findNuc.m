function data = findNuc(data,graphBool,redoBool)
    
    if redoBool == 0 && ~isempty(data.nucPoints)
        return
    end
    
    if isempty(data.nlocs)
        data = chooseSections(data,graphBool,1);
    end
    
    pTime = data.plocsAll{1};
    cutTop=median(pTime);
    pTime(pTime>1.2*cutTop)=cutTop;
    sampleSize = 2 + 3;
    if 0%strfind(data.name,'2-10-12')
        nucleus = 24:38;%10:24;
        data.nucleus = nucleus;
        data.nucSize = 15;
        cytosol = 40:(size(data.xrng,2)-2);
    elseif strfind(data.name,'3-10-12-11')
        nucleus=83:135;
        data.nucleus = nucleus;
        data.nucSize = size(nucleus,2);
        data.cytosol = 140:(size(data.xrng,2)-2);
    elseif strfind(data.name,'2-10-12-2')
        nucleus=99:144;
        data.nucleus = nucleus;
        data.nucSize = size(nucleus,2);
        data.cytosol = 2:95;
    elseif strfind(data.name,'2-10-12-3')
        nucleus=64:140;
        nucleus2=245:331;
        data.nucleus = nucleus;
        data.nucleus2 = nucleus2;
        data.nucSize = size(nucleus,2);
        data.nucSize2 = size(nucleus2,2);
        data.cytosol = 145:240;    
    elseif strfind(data.name,'2-10-12-6')
        nucleus=266:345;
        nucleus2=57:125;
        data.nucleus = nucleus;
        data.nucleus2 = nucleus2;
        data.nucSize = size(nucleus,2);
        data.nucSize2 = size(nucleus2,2);
        data.cytosol = 130:260;    
    elseif strfind(data.name,'2-10-12-7')
        nucleus=70:109;
        data.nucleus = nucleus;
        data.nucSize = size(nucleus,2);
        data.cytosol = 114:212;
    elseif strfind(data.name,'3-10-12-10')
        nucleus=80:142;
        data.nucleus = nucleus;
        data.nucSize = size(nucleus,2);
        data.cytosol = 145:266;
    elseif strfind(data.name,'3-10-12-11')
        nucleus=76:134;
        data.nucleus = nucleus;
        data.nucSize = size(nucleus,2);
        data.cytosol = 140:430;
    elseif strfind(data.name,'3-10-12-1')
        nucleus=91:132;
        data.nucleus = nucleus;
        data.nucSize = size(nucleus,2);
        data.cytosol = 5:86;
    
    elseif 0%strfind(data.name,'28-9-5')
        nucleus = 74:95;
        nucleus2 = 25:42;
        data.nucleus = nucleus;
        data.nucSize = 15;
        nucAvoid = [(nucleus(1)-sampleSize):(nucleus(end)+sampleSize) (nucleus2(1)-sampleSize):(nucleus2(end)+sampleSize)];
        cytosol = setdiff(5:(size(data.xrng,2)-5),nucAvoid);
        
        if graphBool ~= 0
        figure
        subplot(2,1,1)
        plot(data.xrng*data.spaceStep,data.timeStep*pTime)
        xlabel('space (\mu m)')
        ylabel('time (ms)')
        title(strcat('Locations to measure nuclear calcium in',{' '},data.name))
        hold on
        % Plot nucleus in blue
        plot(data.xrng(nucleus)*data.spaceStep, ...
            pTime(nucleus)*data.timeStep,'ob')
        plot(data.xrng(cytosol)*data.spaceStep,pTime(cytosol)*data.timeStep,'oc')
        plot(data.xrng(nucleus2)*data.spaceStep,pTime(nucleus2)*data.timeStep,'om')
        end
        
    else
%     nucThreshold = 50;
    data.nucSize = 0;
    attempt = 0;
%     nucPeak = min([data.nucPeak,pTime]);
    nucPeak=cutTop;
    while data.nucSize<10
        aboveThreshold  = [0 pTime>(nucPeak/2+min(pTime)/2-attempt*10) 0];
        assignin('base','aboveThresh',aboveThreshold);
        % Old version
        % aboveThreshold  = [0 data.plocsAll{1}>(data.nucPeak-nucThreshold/data.timeStep) 0];
        edges           = diff(aboveThreshold);
        rising          = find(edges==1);
        falling         = find(edges==-1);
        nucBegin        = find(rising <= data.nlocs);
        if isempty(nucBegin)
            nucBegin = 10;
        else
            nucBegin = rising(nucBegin(end));
        end
        nucEnd          = find(falling > data.nlocs-1);
        if isempty(nucEnd)
            nucEnd = (size(data.plocsAll{1},2)-10);
        else
            nucEnd = falling(nucEnd(1))-1;
        end

        nucleus = nucBegin:nucEnd;
        data.nucSize = max(size(nucleus));
        attempt = attempt+1;
        if attempt>15
            warning('Not enough points in nucleus');
            figure
            plot(pTime)
            hold on
            refline([0 (nucPeak/2+min(pTime)/2-attempt*10)])
            title(strcat('Not enough points in nucleus:',{' '},data.name))
            break
        end
    end
    figure
    plot(pTime)
    hold on
    refline([0 (nucPeak/2+min(pTime)/2-attempt*10)])
    title('Between attempt switch')
    
    
        attempt = attempt - 2;
    while (data.nucSize>54 && ~isempty(data.nucSize))
        
        aboveThreshold  = [0 pTime>(nucPeak/2+min(pTime)/2-attempt*10) 0];
        assignin('base','aboveThresh',aboveThreshold);
        % Old version
        % aboveThreshold  = [0 data.plocsAll{1}>(data.nucPeak-nucThreshold/data.timeStep) 0];
        edges           = diff(aboveThreshold);
        rising          = find(edges==1);
        falling         = find(edges==-1);
        nucBegin        = find(rising < data.nlocs);
        if isempty(nucBegin)
            nucBegin = 10;
        else
            nucBegin = rising(nucBegin(end));
        end
        nucEnd          = find(falling > data.nlocs-1);
        if isempty(nucEnd)
            nucEnd = (size(data.plocsAll{1},2)-10);
        else
            nucEnd = falling(nucEnd(1))-1;
        end

        nucleus = nucBegin:nucEnd;
        data.nucSize = max(size(nucleus));
        
        figure
        plot(data.xrng*data.spaceStep,data.timeStep*pTime, ...
           data.xrng(nucleus)*data.spaceStep, ...
            pTime(nucleus)*data.timeStep,'ob')
        hold on
        plot([data.xrng(1)*data.spaceStep data.xrng(end)*data.spaceStep], ...
            data.timeStep*[(nucPeak/2+min(pTime)/2-attempt*10) (nucPeak/2+min(pTime)/2-attempt*10)],'r')
        xlabel('space (\mu m)')
        ylabel('time (ms)')
        title(strcat('Locations to measure nuclear calcium in',{' '},data.name))
        
        attempt = attempt-1;
        
        if attempt<-15
            warning('Too many possible points in nucleus');
            figure
            plot(pTime)
            hold on
            refline([0 (nucPeak/2+min(pTime)/2-attempt*10)])
            title(strcat('Too many possible points in nucleus:',{' '},data.name))
            break
        end
    end

    % Double check nucleus
%     crossThreshold = 2*max(data.smoothedImg(data.xrng,nucPeak+500))/3+1/3;
    nucleus2=nucleus;
    offset = 1000;
    line = data.smoothedImg(data.xrng,data.init+offset);
    crossThreshold = sort(line,'ascend');
    crossThreshold = crossThreshold(round(0.65*size(crossThreshold,1)));
    nCross = find(line>crossThreshold);
    nucleus = intersect(nCross,nucleus);
    nucInt = [0 cumsum(diff(nucleus')~=1)];
    nucleus = nucleus(nucInt==mode(nucInt));
    
    data.nucleus = nucleus;
    data.nucSize = max(size(data.nucleus));
%     figure
%     plot(line)
%     hold on
%     refline([0 crossThreshold])
%     plot(data.nucleus,line(data.nucleus),'ob')
%     for jj=1:10
%         line = data.smoothedImg(data.xrng,data.init+offset+20*jj);
%         crossThreshold = sort(line,'ascend');
%         crossThreshold = crossThreshold(round(0.7*size(crossThreshold,1)));
%         plot(line)
%         refline([0 crossThreshold])
%     end
    
    
    
%     m = [4,4,4,10,10,10,4,4,5];
% 
% [U,~,i1]=unique(m);
% freq= histc(m,U);
% [~,i2] = sort(freq(i1),'descend');
% 
% m(i2)
%     falling                 = find(edges==1);
%     rising               = find(edges==-1);
    upSizes             = falling-rising;
    [~, maxULoc]        = max(upSizes);
    maxU2Loc            = maxULoc;
%     [~, maxU2Loc]       = max([upSizes(1:maxULoc-1) upSizes(maxULoc+1:end)]);
%     if maxU2Loc >= maxULoc
%         maxU2Loc = maxU2Loc+1;
%     end
%     maxULoc
%     maxU2Loc
    cytSections = {sampleSize:min([rising(maxULoc),rising(maxU2Loc)])-1, ...
        min([falling(maxULoc),falling(maxU2Loc)]):max([rising(maxULoc),rising(maxU2Loc)])-1, ...
        max([falling(maxULoc),falling(maxU2Loc)]):(size(data.xrng,2)-sampleSize)};
    [~,maxLoc] = max([min([rising(maxULoc),rising(maxU2Loc)])-1; ...
        max([rising(maxULoc),rising(maxU2Loc)])-min([falling(maxULoc),falling(maxU2Loc)]); ...
        size(aboveThreshold,2)-max([falling(maxULoc),falling(maxU2Loc)])]);
    cytosol = cytSections{maxLoc};
    cytosol = cytosol(5):cytosol(end-5);
     if graphBool ~= 0
        figure
        subplot(2,1,1)
        plot(data.xrng*data.spaceStep,data.timeStep*pTime)
%         plot(data.xrng(data.nlocs)*data.spaceStep,data.nucPeak*data.timeStep,'or')
        xlabel('space (\mu m)')
        ylabel('time (ms)')
        title(strcat('Locations to measure nuclear calcium in',{' '},data.name))
        hold on
%         plot(data.xrng(data.clocs)*data.spaceStep,cyt*data.timeStep,'og')
        % Plot nucleus in blue
        plot(data.xrng(nucleus)*data.spaceStep, ...
            pTime(nucleus)*data.timeStep,'ob')
        plot(data.xrng(cytosol)*data.spaceStep,pTime(cytosol)*data.timeStep,'oc')
        % Plot 2nd biggest area in green
        plot(data.xrng(rising(maxU2Loc):falling(maxU2Loc)-1)*data.spaceStep, ...
            pTime(rising(maxU2Loc):falling(maxU2Loc)-1)*data.timeStep+10,'og')
        plot(data.xrng(rising(maxULoc):falling(maxULoc))*data.spaceStep,pTime(rising(maxULoc):falling(maxULoc))*data.timeStep+10,'om')
        plot(data.xrng(nucleus2)*data.spaceStep, ...
            pTime(nucleus2)*data.timeStep,'or')
        plot(data.xrng(nCross)*data.spaceStep, ...
            pTime(nCross)*data.timeStep,'oy')
        hold off
        subplot(2,1,2)
        plot(data.xrng*data.spaceStep,aboveThreshold(2:end-1))
        if data.saveBool ~=0
            saveas(gcf,strcat(data.fpath,'graphs/cytnuc-',data.name,'.png'))
        end
    end
    
    end
%     if upSizes((maxULoc+1):
    
% 
%     cytBegin        = find(rising>data.nlocs);
%     if isempty(cytBegin)
%         cytosol = falling(nucEnd(1))+2:(size(data.plocsAll{1},2)-10);
%     else
%         cytosol = falling(nucEnd(1))+2:rising(cytBegin(1))-2;
%     end
%     
%     if nucEnd(1)>1
%         cytosol2 = falling(nucEnd(1)-1)+2:rising(nucBegin(end))-2;
%     else
%         cytosol2 = 10:rising(nucBegin(end))-2;
%     end
%     
%     if size(cytosol2,2)>size(cytosol,2)
%         cytosol = cytosol2;
%     end
    

    
   
    
    if(size(data.cytosol,2)<9) || max(size(data.nucleus))<5
        warning(strcat('Cytosol:',num2str(size(data.cytosol,2)),' or nucleus:',num2str(max(size(data.nucleus))),' too small'));
        return
    end
    % Choose sections of nucleus and cytosol
    data.nucPoints = datasample(data.nucleus,5,'Replace',false);
    data.cytPoints = datasample(data.cytosol(3:end-2),5,'Replace',false);
%     data.pnrPoints = [data.nucleus(1)-2, data.nucleus(end)+2];
end
