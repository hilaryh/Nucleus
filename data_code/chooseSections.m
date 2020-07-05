function data = chooseSections(data,graphBool,redoBool)    
    if redoBool == 0 && ~isempty(data.clocs)
        return
    end
    if isempty(data.smoothedImg)
       data = processImg(data,0,1);
    end
    if isempty(data.xrng)
       data = findScanBorder(data,0,1);
    end
    if isempty(data.plocsAll)
        data = findCriticalPoints(data);
    end

    nucThreshold = 10;
    
    [nucs,nlocs] = findpeaks(data.plocsAll{1}(nucThreshold:(end ...
        -nucThreshold)),'MinPeakDistance',25);
    nlocs        = nlocs + nucThreshold-1;
    [data.nucPeak,maxPos]= max(nucs);
    data.nlocs        = nlocs(maxPos);

    [cyt,clocs] = findpeaks(-data.plocsAll{1}(nucThreshold:( ...
        size(data.xrng,2)-nucThreshold)),'MinPeakDistance',25);
    clocs       = clocs + nucThreshold - 1;
    cyt         = -cyt;
    [cyt,minPos]= min(cyt);
    data.clocs       = clocs(minPos);

%     [gmax,gtime] = max(gradient(data.smoothedImg(data.clocs+data.xrng(1)-2,data.init:(data.init+2500))));
%     data.init = gtime+data.init-19
%     figure
%     plot((data.init-100):(data.init+2500), ...
%         data.smoothedImg(data.clocs+data.xrng(1)-2,(data.init-100):(data.init+2500)), ...
%         data.init,data.smoothedImg(data.clocs+data.xrng(1)-2,data.init),'or')
    
    if graphBool == 1
        figure
        plot(data.xrng*data.spaceStep,data.timeStep*data.plocsAll{1}, ...
            data.xrng(data.nlocs)*data.spaceStep,data.nucPeak*data.timeStep,'or')
        xlabel('space (\mu m)')
        ylabel('time (ms)')
        title(strcat('Locations to measure nuclear calcium in',{' '},data.name))
        hold on
        plot(data.xrng(data.clocs)*data.spaceStep,cyt*data.timeStep,'og')
        plot([data.xrng(nucThreshold) data.xrng(nucThreshold)]*data.spaceStep,...
            [min(data.plocsAll{1}) max(data.plocsAll{1})]*data.timeStep)
        plot([data.xrng(size(data.xrng,2)-nucThreshold) ...
            data.xrng(size(data.xrng,2)-nucThreshold)]*data.spaceStep, ...
            [min(data.plocsAll{1}) max(data.plocsAll{1})]*data.timeStep)
        hold off
        if data.saveBool == 1
            saveas(gcf,strcat(data.fpath,'graphs/release-',data.name,'.png'))
        end
    end
end
