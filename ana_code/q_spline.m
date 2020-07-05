% H Hunt 2019
% Generate boundary condition q to simulate cytosolic calcium in a
% cardiomyocyte with provided parameters.
% t given in seconds, q calculated in uM
% Called by avCyt.m in cutNuc1DData.m
function sol=q_spline(t,fits,crtpts)
% global peakMin;global peakMax
if size(fits,1)>1
    fits=fits(1,:);
    crtpts=crtpts(1,:);
end

% Work out which fit applies
t1=t(t<=fits{1}.breaks(1));
t2=t(t>fits{1}.breaks(1)&t<=fits{1}.breaks(end));
t3=t(t>fits{1}.breaks(end));
% sol1=fits{2}(crtpts(end))*ones(size(t1));
sol1=ppval(fits{1},fits{1}.breaks(1))*ones(size(t1));
sol2=ppval(fits{1},t2);
sol3=ppval(fits{1},fits{1}.breaks(end))*ones(size(t3));

sol=[sol1,sol2,sol3];
% sol=(sol*(peakMax-peakMin)-0.1*peakMax+peakMin)/0.9;
end