function dilateFit = dilatepp(inputFit,inclength,incheight)
newBreaks=[inputFit.breaks(1:7),(inputFit.breaks(8:end)-inputFit.breaks(7))*inclength+inputFit.breaks(7)];
newVals=incheight*(ppval(inputFit,inputFit.breaks)-ppval(inputFit,inputFit.breaks(1)))+ppval(inputFit,inputFit.breaks(1));
f = fit(newBreaks',newVals','smoothingspline','SmoothingParam',0.07);
dilateFit=coeffvalues(f);
end