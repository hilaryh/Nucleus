function [dataStruct] = loadData(parentFolder,nucs,dataSheet)
    dataStruct = dir(strcat(parentFolder,'*.tif'));
    dataStruct = rmfield(dataStruct,{'isdir','datenum','bytes'});
    dataStruct = dataStruct(nucs);
    dataStruct(1).Dvals=[];
    
    for ii=1:size(dataStruct,1)
        dataStruct(ii).timeStep   = dataSheet{nucs(ii),16}; %ms
        dataStruct(ii).spaceStep  = dataSheet{nucs(ii),14}; %um
        dataStruct(ii).imgO = imread(strcat(parentFolder,dataStruct(ii).name));
        dataStruct(ii).imgO = double(dataStruct(ii).imgO);
        dataStruct(ii).name       = dataStruct(ii).name(1:end-4);
    end
end