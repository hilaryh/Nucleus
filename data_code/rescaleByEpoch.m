function newImg = rescaleByEpoch(data)
cytP = data.cytPoints(1);
rawline = sum(data.img((cytP-2):(cytP+2),:))/5;
[line, ~] = bior11tot1(rawline);
line=line(:,data.trng);
[rpks, rlocs] = findpeaks(line,'MinPeakProminence',0.2,'MinPeakDistance',2000);
newImg = data.img;
types = find(diff(rlocs)>6000);
assignin('base','atypes2',rpks)
assignin('base','alocs2',rlocs)
if size(types,2)>3
    newd = diff([0 types]);
    rpks(types(newd==1))=[];
    rlocs(types(newd==1))=[];
    types(newd==1)=[];
    while sum(newd==2)>0

        xtypes = find(newd==2);
        types(xtypes(1))=[];
        newd = diff([0 types]);
    end
    while types(end)==size(rlocs,2)
        types(end)=[];
        rlocs(end)=[];
        rpks(end)=[];
    end
end
if size(types,2)<3
    types(3) = size(rlocs,2);
end
if rpks(types(3))<rpks(types(3)-1)-0.3
    types(3)=types(3)-1;
end

[bpks, blocs] = findpeaks(-line,'MinPeakProminence',0,'MinPeakDistance',2000);
bpks=-bpks;

%second epoch
for ii = 1:2
f0 = find(blocs<rlocs(types(ii)+1),1,'last');
f0=blocs(f0);
mean(mean(newImg(:,(f0-5):f0)))
newImg(:,(f0-5):f0) = mean(mean(newImg(:,(f0-5):f0)))*ones(size(newImg(:,(f0-5):f0)));
newImg(:,f0:end) = fonf0(newImg(:,f0:end));
end
newline = newImg(data.cytPoints(1),data.trng);
rawline = data.img(data.cytPoints(1),data.trng);

% figure
% plot(line)
% hold on
% plot(rlocs, rpks,'oy')
% plot(rlocs(types),3*ones(size(types)),'og')
% plot(blocs,bpks,'ob')
% plot(f0,1,'xr')
% 
% figure
% plot(newline)
% hold on
% plot(rawline)


end