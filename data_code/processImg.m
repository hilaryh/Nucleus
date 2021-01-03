function data = processImg(data,graphBool,redoBool)
    if redoBool == 0 && ~isempty(data.smoothedImg)
        return
    end
    
    % Using automatically generated denoising functions 
    % from the SWT Denoising 1D toolbox
	maxTime=15000
    img = data.imgO;
    place = 200*size(img,1)/512;
%     minProm = 0.3*(max(max(img))-min(min(img)));
    [~, ~,~,proms] = findpeaks(mean(img(:,1:maxTime)));
    minProm=max(proms)/2;
    [~, plocs] = findpeaks(mean(img(:,1:maxTime)),'MinPeakProminence',minProm);
%     figure
%     plot(img(place,1:3500))
%     hold on 
%     plot(plocs,pks,'or')
%     hold off
%     figure
%     subplot(2,1,1)
%     plot(img(place,:))
    img = flattenInit(img,plocs(1)-40);
    data.initFluo = mean(img(1,:));
    
    img=fonf0(img);
%     subplot(2,1,2)
%     plot(img(place,:))
    newImg = zeros(size(img));
    
    if size(newImg,2)> 4e5
        parfor ii = 1:size(img,1)
             [smoothLine, ~] = bior11tot(img(ii,:));
             maxCa=10*median(smoothLine);
             ex = find(smoothLine>maxCa | smoothLine<0);
             tot = [ex ex+1 ex-1];
             smoothLine(ex)=median(smoothLine);
             newImg(ii,:)=smoothLine;
        end
    else
                parfor ii = 1:size(img,1)
             [smoothLine, ~] = bior11tot1(img(ii,:));
             maxCa=10*median(smoothLine);
             smoothLine(smoothLine>maxCa | smoothLine<0)=median(smoothLine);
             newImg(ii,:)=smoothLine;
        end
    end
        
	data.smoothedImg = newImg;
    data.img=img;
    newImg = data.smoothedImg;
    
    if graphBool==1
        figure
        plot(img(place,:),'Color',[1 0.4 0.3])
        hold on
        plot(newImg(place,:),'b')
    end
    
end
