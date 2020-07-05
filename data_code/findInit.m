    % Set init to just before the calcium transient appears
    
    function data = findInit(data,graphBool,redoBool)
    if redoBool == 0 && ~isempty(data.init)
        return
    end
    newImg=data.smoothedImg;
    img = data.img;
    spacesize=size(img,1);
    xrng=round(2*spacesize/5:4*spacesize/5);
    peakSeparation = 1500;
    fpLine=median(newImg(xrng,1:15000));
    assignin('base','fpLine',fpLine)
    [~, ~,~,proms] = findpeaks(fpLine);
    minProm=max(proms)/2;
    data.minProm=minProm;
    [apks, alocs] = findpeaks(fpLine,'MinPeakDistance',peakSeparation,'MinPeakProminence',minProm);
    [pks, plocs] = findpeaks(-fpLine);
    if size(alocs,2)>=2
        if apks(2)>1.2*apks(1)
            apks(1)=[];
            alocs(1)=[];
        end
    end
    if pks(2)>pks(1)+0.05
        pks(1)=[];
        plocs(1)=[];
    end
    alocs2 = find(plocs<alocs(1));
    if isempty(alocs2)
        alocs(1)=[];
        apks(1)=[];
        alocs2 = find(plocs<alocs(1));
    else
        alocs2(pks(alocs2)<-apks+0.05)=[];
    end

    data.plocs = plocs(alocs2(end));
    pks = pks(alocs2(end));
%     figure
%     plot(-newImg(place,1:3500))
%     hold on
%     plot(data.plocs,pks,'or')
%     plot(alocs,-apks*0.9,'ob')
    data.init = max(round(data.plocs(1)),alocs(1)-700);
    if data.init<0
        data.init = 10;
        warning('Not enough time before peak');
    end
    
    % Rewrite the first timestep so that it is the average of the pre-transient
    % fluorescence
    smoothedImg = flattenInit(newImg,data.init);
    img = flattenInit(img,data.init);

    % Rescale the graphs to show f/f0
    data.smoothedImg = fonf0(smoothedImg);
    data.smoothedImg(:,1:20) = ones(size(smoothedImg(:,1:20)));
    data.img = fonf0(img);

    % Account for brightening of fluorescent dye 
    %...
    
    if graphBool==1
        figure
        plot((1:size(smoothedImg,2))*data.timeStep,median(smoothedImg(xrng,:)),data.timeStep* ...
            (data.plocs),-pks,'or')
        xlabel('time (ms)')
        ylabel('calcium (f/f_0)')
        title(strcat('Looking for init at',{' '},data.name))
        hold on
        plot([data.timeStep*(data.init) data.timeStep*(data.init+150) ...
            data.timeStep*(data.init+300)], [median(smoothedImg(xrng,data.init))+1 ...
            median(smoothedImg(xrng,data.init+150)) median(smoothedImg(xrng,data.init+300))],'og')
        hold off
        if data.saveBool == 1
            saveas(gcf,strcat(data.fpath,'graphs/initFind-',data.name,'.png'))
        end
    end
end
