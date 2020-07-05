function data = generateRows(data,row,redoBool)
    if isfield(data,'nucSize2') 
        if ~isempty(data(row).img) && redoBool==0
            return
        end
    end
    data(row).initFluo = [];
%     data(row).imgO = imread(strcat(data(row).fpath,data(row).name,'.tif'));
%     data(row).imgO = double(data(row).imgO);
    data(row).img = [];
    data(row).minProm=[];
    data(row).spaceStep = data(row).spaceStep ...
        *512/size(data(row).imgO,1);
    data(row).init = [];
    data(row).smoothedImg = [];
    data(row).butterworth = [];
    data(row).xrng = [];
    data(row).trng = [];
    data(row).plocsAll = [];
    data(row).plocs = [];
    data(row).pks = [];
    data(row).tlocs = [];
    data(row).clocs = [];
    data(row).nlocs = [];
    data(row).nucPeak = [];
    data(row).nucSize = [];
    data(row).nucSize2 = [];
    data(row).nucleus = [];
    data(row).nucleus2 = [];
    data(row).cytosol = [];
    data(row).nucPoints = [];
    data(row).cytPoints = [];
%     data(row).pnrPoints = [];
%     data(row).ttPeakCyt = [];
%     data(row).ttPeakNuc = [];
%     data(row).ttBaseCyt = [];
%     data(row).ttBaseNuc = [];
%     data(row).maxCyt = [];
%     data(row).maxNuc = [];
%     data(row).maxCytB = [];
%     data(row).maxNucB = [];
%     data(row).cfwhm = [];
%     data(row).cfwhmB = [];
%     data(row).nfwhm = [];
%     data(row).calcLevel=[]; 
%     data(row).linesC=[]; 
end