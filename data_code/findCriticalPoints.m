% Look for peaks and troughs in data
function data = findCriticalPoints(data,redoBool)
    if redoBool == 0 && ~isempty(data.plocsAll)
        return
    end
    if isempty(data.smoothedImg)
       data = processImg(data,0,redoBool);
    end
    if isempty(data.xrng)
       data = findScanBorder(data,0,redoBool);
    end
    
    peakSeparation = 1500;
    pks = cell(size(data.xrng,2),1);
    % The time of calcium peaks at every point
    data.plocs = cell(size(data.xrng,2),1);
    for ii=1:size(data.xrng,2)
        [pks{ii}, data.plocs{ii}] = ...
            findpeaks(data.smoothedImg(ii+data.xrng(1)-1,data.init:end), ...
            'MinPeakDistance',peakSeparation,'MinPeakProminence',data.minProm);
        data.plocs{ii} = data.plocs{ii}+data.init-1;
%         if mod(ii,10)==0
%         figure
%         plot(data.smoothedImg(ii+data.xrng(1)-1,:), ...
%             (data.plocs{ii})*data.timeStep,pks{ii},'or')
%         end
    end
    % The peak calcium levels at every point
    data.pks = pks;
    % Find troughs
    trs = cell(size(data.xrng,2),1);
    data.tlocs = cell(size(data.xrng,2),1);
    for ii=1:size(data.xrng,2)
        [trs{ii},data.tlocs{ii}] = ...
            findpeaks(-data.smoothedImg(ii+data.xrng(1)-1,:), ...
            'MinPeakDistance',peakSeparation);
        trs{ii}=-trs{ii};
        if size(data.tlocs{ii},2)<2
            data.tlocs{ii}(2) = size(data.smoothedImg,2);
            trs{ii} = data.smoothedImg(ii+data.xrng(1)-1,data.tlocs{ii}(2));
        end
    end
    
    % plocsAll is a cell array consisting of a single cell with an array
    % containing the location of the first peak in calcium at each point in
    % the myocyte
    
    % I forget why it's not just an array. The cell part is unneccessary at
    % the moment (24/3/2017).
    data.plocsAll=cell(1,1);
    for jj = 1:size(data.plocsAll,1)
        data.plocsAll{jj}=[];
        for ii=1:size(data.xrng,2)
            if size(data.plocs{ii},2)<jj
                data.plocsAll{jj}=[data.plocsAll{jj} NaN];
            else
                data.plocsAll{jj}=[data.plocsAll{jj} data.plocs{ii}(jj)];
            end
        end
    end
    if max(data.plocsAll{1})>1.8*min(data.plocsAll{1})
        firstPts=data.plocsAll{1};
        oddPts = find(firstPts>1.8*min(firstPts));
        oddPtsLo = oddPts(oddPts<size(data.xrng,2)/3);
        if isempty(oddPtsLo)
            oddPtsLo=0;
        end
        oddPtsHi = oddPts(oddPts>2*size(data.xrng,2)/3);
        if isempty(oddPtsHi)
            oddPtsHi=size(data.xrng)+1;
        end
        if isequal(oddPts,union(oddPtsHi,oddPtsLo))
            data.xrng=data.xrng((max(oddPtsLo)+1):(min(oddPtsHi)-1));
            data.plocsAll{1}=data.plocsAll{1}((max(oddPtsLo)+1):(min(oddPtsHi)-1));
        else
            warning('Some parts of the cell are peaking far too late')
        end
    end
            
        
    figure
    plot(data.plocsAll{1});
    title('findCritPts: plocsAll')
end
