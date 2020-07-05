function data = butterworthImg(xdata,redoBool)
data = xdata;
if redoBool == 0 && ~isempty(xdata.butterworth)
    return 
end
    data.butterworth = filterNoisyLines(data.img,1000);
    
    if redoBool==1
        figure
        plot(data.img(50,:),'Color',[0.8 1 1])
        hold on
        plot(data.butterworth(50,:),'Color','b')
        plot(data.smoothedImg(50,:),'--','Color','g')
    end
end
