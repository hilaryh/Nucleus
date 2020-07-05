% function cutNuc1DData(dataStruct)
onserver=0;
if exist('/home/hhunt1/AnaNuc/1Dsims/data/','dir')
    onserver=1;
end
% On laptop
if 1-onserver
    homeDir='/Volumes/UniWD/Uni/data_Bass/';
    processDir='/Users/hhunt1/Documents/Write Ups/Maths/Data/1D/';
    addpath('/Users/hhunt1/Documents/Write Ups/Maths/Data/cutNuc1DData/')
else % On server
    homeDir='/home/hhunt1/AnaNuc/1Dsims/';
    processDir='/home/hhunt1/AnaNuc/1Dsims/data/';
    addpath('/home/hhunt1/AnaNuc/1Dsims/numerical/')
end
%% Load info from spreadsheet
    [~,~,dataSheet] = xlsread(strcat(homeDir,'linescan_data_guide.xlsx'));
    nucCol=[dataSheet{2:end,10}];
    nucs=find(nucCol==1)+1;
    nucAngles=NaN*ones(size(nucs,2),2);
    anglesCol(:,1)=dataSheet(nucs,9);
    dblnuc=cellfun(@(x)find(x=='_'),anglesCol,'Un',0);
    for ii=1:size(dblnuc,1)
        if dblnuc{ii}
            nucAngles(ii,2)=str2double(anglesCol{ii}((dblnuc{ii}+1:end)));
            nucAngles(ii,1)=str2double(anglesCol{ii}(1:(dblnuc{ii}-1)));
        else
            nucAngles(ii,1)=anglesCol{ii};
        end
    end
% Identify which cells have nucleus info
    cellID=dataSheet(nucs,3);
%%    Analyse cells
    redoBool = 0;
    saveBool=0;
    graphBool=ones(7,1);
    prepath=strcat(processDir,'data','.mat');
    if isempty(dataStruct)
        if exist(prepath,'file')==2
            load(prepath)
        else
            [dataStruct]=loadData(strcat(homeDir,'data_montage/'),nucs,dataSheet);
%                 dataStruct=identifyNuc(dataStruct);
%                 save(prepath,'dataStruct');
        end
    end
    for cellit=2:size(cellID,1)
        
        dataStruct(cellit).fpath = strcat(homeDir,'data_montage/');
        dataStruct = generateRows(dataStruct,cellit,redoBool);
        dataStruct(cellit).fpath = processDir;
        dataStruct(cellit).saveBool = saveBool;
        dataStruct(cellit) = processImg(dataStruct(cellit),graphBool(1),0);
        dataStruct(cellit) = findInit(dataStruct(cellit),(1),0);
    %     dataStructTt(cellit) = processImg(dataStructTt(cellit),graphBool(1),redoBool);
        dataStruct(cellit) = butterworthImg(dataStruct(cellit),0);
        % Find borders of cell to crop linescan data
        % Fill in fields for xrng, trng
        dataStruct(cellit) = findScanBorder(dataStruct(cellit),graphBool(2),redoBool);
        % Look for peaks and troughs in data
        % Fill in fields plocs,pks,tlocs,plocsAll
        % Apart from pks, these are looking at timing of the calcium transient
        dataStruct(cellit) = findCriticalPoints(dataStruct(cellit),redoBool);

        % Identify sections of nucleus and cytosol
        dataStruct(cellit) = chooseSections(dataStruct(cellit),graphBool(3),redoBool);
        % Identify nucleus
        % Fill in fields nucSize,nucleus
        dataStruct(cellit).cytosol=[];
        dataStruct(cellit) = findNuc(dataStruct(cellit),graphBool(4),1);
%         dataStruct(cellit).smoothedImg = rescaleByEpoch(dataStruct(cellit));
        topSurfGraph(dataStruct(cellit),1,2)
%         save(prepath,'dataStruct','-v7.3');
        cutangle=dataSheet(nucs(cellit),11);
        cutangle=cutangle{1};
        cutangle=cutangle*pi/180;
        dataStruct(cellit).cutangle=cutangle;
        
        [dataStruct(cellit).crtpts,dataStruct(cellit).mccpos,dataStruct(cellit).fstfit,dataStruct(cellit).thisPeak]=...
            avCyt_spline(dataStruct,cellit,1);
        
        extra=10*dataStruct(cellit).spaceStep;
        cutLabel=mod(dataStruct(cellit).cutangle*180/pi,180);
        if cutLabel>90
            cutLabel=180-cutLabel;
        end
        mtpath=strcat(processDir,'mt',num2str(cutLabel),...
            '-',num2str(4.7),'-',num2str(18),'.mat');
%         save(prepath,'dataStruct','-v7.3');
        if exist(mtpath,'file')==2
            load(mtpath)
        else
            [dataStruct(cellit).maxtime,dataStruct(cellit).fdhmtime,dataStruct(cellit).Drange,dataStruct(cellit).mtimes]=...
                findDCurve_spline_num(dataStruct,cellit,2,1);
            maxtime=dataStruct(cellit).maxtime;
            save(mtpath,'maxtime','-v7.3');
        end
        [dataStruct(cellit).Dvalsm,dataStruct(cellit).Dvalsfdhm,dataStruct(cellit).fdhm,dataStruct(cellit).ttpeak,dataStruct(cellit).nfdhm,dataStruct(cellit).nttpeak]...
            =fitDinNuc_spline_num(dataStruct,cellit,2,1);
        save(prepath,'dataStruct','-v7.3');
    end