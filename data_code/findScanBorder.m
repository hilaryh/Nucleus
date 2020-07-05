function data = findScanBorder(xdata,graphBool,redoBool)
    data = xdata;
    if redoBool == 0 && ~isempty(xdata.xrng)
        return
    end

    if isempty(data.img)
       data = dataInit(data,0,1);
    end   
    timePoint = data.init+400;
    [pks, locs] = ...
        findpeaks(data.smoothedImg(:,timePoint));
    
    % Check for absurdy bright sections
    [spks,slocs]=(sort(pks));
    pkdiffs=diff(spks);
    numpks=size(pkdiffs,1);
    oddity=find(pkdiffs(floor(numpks/2):end)>2)+floor(numpks/2);
    pks(slocs(oddity:end))=[];
    locs(slocs(oddity:end))=[];
    
    cutOff = 1.05;
    if max(pks)<cutOff
        max(pks)
        figure 
        plot(data.smoothedImg(1:end,timePoint))
        return
    end
    ii=1;
    while pks(ii)<cutOff
        ii=ii+1;
    end
%     pks(ii)
    xmin = locs(ii)-1;
    
    jj=size(pks,1);
    while pks(jj)<cutOff
        jj=jj-1;
    end
    xmax = locs(jj)-1;
    
    tmin = max(data.init-600,1);      % the beginning
    tmax = size(data.smoothedImg,2);  
    
    data.xrng = (xmin:xmax);
    data.trng = tmin:tmax;
    if graphBool == 1
        figure
        xpos = (0:size(data.smoothedImg,1)-1)*data.spaceStep;
        plot(xpos,data.smoothedImg(:,timePoint),xpos(locs),pks,'or')
        title(strcat('Finding edge of cell in',{' '},data.name,' at', ...
            num2str((timePoint)*data.timeStep),{' '},num2str(ii),{' '},num2str(jj)))
        xlabel('space (\mu m)')
        ylabel('fluorescence (F/F_0)')
        hold on
        plot([0 xmax*data.spaceStep],[cutOff cutOff])
        plot([xmin xmin]*data.spaceStep, [min(data.smoothedImg(:,timePoint)) max(data.smoothedImg(:,timePoint))])
        plot([xmax xmax]*data.spaceStep, [min(data.smoothedImg(:,timePoint)) max(data.smoothedImg(:,timePoint))])
        hold off
        if data.saveBool == 1
            saveas(gcf,strcat(data.fpath,'graphs/edgeFind-',data.name,'.png'))
        end
    end
end
