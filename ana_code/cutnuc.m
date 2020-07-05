function [sol,h]=cutnuc(cutangle,l,a,sol2d,t)
opponadj=tan(cutangle);
if opponadj<0
    cutangle=pi-mod(cutangle,pi);
    opponadj=tan(cutangle);
end
cutlength=201;
interval=(1/cutlength);
newa=size(sol2d{1},1);
newl=floor(size(sol2d{1},2)/2);
if opponadj==0
    hyplen=l;
    rh=ones(1,cutlength+1);
    zh=round(newl*(interval:interval:1))+newl;
elseif a/l<=opponadj
%     warning('redo this section of if statement')
    hyplen=sqrt(a^2+(a/opponadj)^2);
%     newa=a;
%     newl=a/opponadj;
    rh=round(newa*(interval:interval:1));
    zh=round(newl*a*(interval:interval:1)/(l*opponadj))+newl;
else
    hyplen=sqrt((l)^2+(l*opponadj)^2);
%     newa=l*opponadj;
%     newl=l;
    rh=round(newa*l*opponadj*(interval:interval:1)/a);
    zh=round(newl*(interval:interval:1))+newl;
end
rh(rh==0)=1;
zh(zh==0)=1;
h=0:interval:hyplen;
sol=cell(size(t));
if ~isempty(zh)
    for time=1:size(t,2)
%         linepic=sol2d{50};
        try
%         sol{time}=tsol{time}(idx([end:-1:2 2 2:end]),idy)+integral(fun,0,t(time));
for hstep=1:cutlength
        sol{time}(hstep,1)=sol2d{time}(rh(hstep),zh(hstep));%integral(fun,0,t(time));
%         linepic(rh(hstep),zh(hstep))=0;
end
%         sol{time}=rsol(idx([end:-1:2 2 2:end]))...
%             .*zsol(idy)*tsol{time}+integral(fun,0,t(time));
        catch
            warning('oh no! :(')
        end
    end
else
    for time=1:size(t,2)
        warning('add zcoeff at line 90,sim_nuc_cyl')
        sol{time}=sol2d{time}(ones(1,cutlength),ones(1,cutlength));%+integral(fun,0,t(time));
%         sol{time}=rsol(2*ones(1,2*size(h,2)-1))*tsol{time}+integral(fun,0,t(time));
    end
end
end