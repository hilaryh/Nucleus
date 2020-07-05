function [init img minprom place] = findInit1D(dataStructTt)
img=dataStructTt.imgO;
testprom=0.5*(max(img(:))-min(img(:)));
[t1,t2]=findpeaks(mean(img),'MinPeakProminence',testprom);

end