% H Hunt 2019
% Used in findDvals

function [cutca,h,tsol,rfun,zfun,tfun]=sim_nucleus_cyl_wz_spline_cstq(D,l,a,cutangle,tint,tmax)
interval=0.01;
r=0:interval:a;
z=0:interval:l;
% maxn=100;
tmin=0;
% tmax=1;
% tint=0.01;
t=tmin:tint:tmax;
% add padding
%t=[tmin*ones(1,crtpts(1)) t(1:(end-crtpts(1)))];
[rfun,zfun,tfun,tsol]=generatecaconc(a,l,D,r,z,t);
[cutca,h]=cutnuc(cutangle,l,a,r,z,tsol,t);
if 1
    makemovie(tmin,tint,D,tsol)
    makemovie1D(tmin,tint,D,cutca,h,cutangle)
end
end
function [rfun,zfun,just_tfun,tsol]=generatecaconc(a,l,D,r,z,t)
if isempty(D)
    D=[];
end
v_o=2;
v_i=0-v_o;
maxt=1e3;
% maxt=maxt+939-mod(maxt,939);
% maxt=maxt+469-mod(maxt,469)+374;
% maxt=20*940;
tsol=cell(size(t));
alphan=zeros(1,maxt+1);
for n=1:maxt
    alphan(n+1) = fzero(@(z) besselj(0, a*z), [alphan(n)+0.1/a, alphan(n)+5/a]);
end
rfun=besselj(0,r'*alphan(2:end))./besselj(1,alphan(2:end)*a)./alphan(2:end);
zfun=sin(z'*(2*(0:(maxt-1))+1)*pi/l)./(2*(0:(maxt-1))+1);
tsol=cellfun(@(x)zeros(size(rfun(:,1)*zfun(:,1)')),tsol,'Un',0);
common_coeff=8/(pi*a);
just_tfun=cell(1,size(t,2));
for time=1:size(t,2)
    n=1:maxt;
    m=0:  (maxt-1);
    if t(time)==0
        just_tfun{time}=common_coeff*v_i*ones(maxt,maxt);
    else
        just_tfun{time}=common_coeff*v_i*exp(-D*(alphan(n+1)'.^2+((2*m+1)*pi/l).^2)*t(time));
    end
    tsol{time}=rfun*just_tfun{time}'*zfun'+v_o;
%     tsol{time}=tsol{time};
    
    if abs(tsol{time})>5
        warning('something is not right with equations')
    end
%     tsol{time}=sum(cat(3,tfun{:}),3);
end
end

function [sol,h]=cutnuc(cutangle,l,a,r,z,tsol,t)
opponadj=tan(cutangle);
if opponadj<0
    cutangle=pi-mod(cutangle,pi);
    opponadj=tan(cutangle);
end
interval=0.01;
if 2*a/l<=opponadj
    hyplen=sqrt(a^2+(a/opponadj)^2);
    h=0:interval:hyplen;
    y=round(0:(a/(size(h,2)-1)):a,2);
    x=round(0:(a/opponadj/(4*(size(h,2)-1))):(a/opponadj/2),2);
    [~,idr] = ismember(y,round(r,2));
    idr=nonzeros(idr);
    [~,idz]=ismember(x,round(z,2));
    idz=nonzeros(idz);
else
    hyplen=sqrt((l/2)^2+(l*opponadj/2)^2);
    h=0:interval:hyplen;
    y=round(0:(l*opponadj/(2*(size(h,2)-1))):(l*opponadj/2),2); % a dir
    x=round(0:(l/(size(h,2)-1)/2):(l),2);
    idz=arrayfun(@(s)find(s==round(z,2),1,'first'),round(x,2),'UniformOutput',false);
    idz=[idz{:}];
    idr=arrayfun(@(s)find(s==round(r,2),1,'first'),round(y,2),'UniformOutput',false);
    idr=[idr{:}];
end
sol=cell(size(tsol));

if ~isempty(y)
    for time=1:size(tsol,2)
        try
%         sol{time}=tsol{time}(idx([end:-1:2 2 2:end]),idy)+integral(fun,0,t(time));
%         Does this need to change if 2*a/l<=opponadj?
        sol{time}=tsol{time}((idz-1)*size(tsol{time},1)+(idr([end:-1:2 2 2:end])));%integral(fun,0,t(time));
%         sol{time}=rsol(idx([end:-1:2 2 2:end]))...
%             .*zsol(idy)*tsol{time}+integral(fun,0,t(time));
        catch
            warning('oh no!')
        end
    end
else
    for time=1:size(tsol,2)
        warning('add zcoeff at line 90,sim_nuc_cyl')
        sol{time}=tsol{time}(2*ones(1,2*size(h,2)-1),ones(1,2*size(h,2)-1))+integral(fun,0,t(time));
%         sol{time}=rsol(2*ones(1,2*size(h,2)-1))*tsol{time}+integral(fun,0,t(time));
    end
end
end

function makemovie1D(tmin,tint,D,cutca,h,cutangle)
if exist('/home/hhunt1/AnaNuc/1Dsims/data/','dir')
    savefolder='/home/hhunt1/AnaNuc/1Dsims/data/';
else
    savefolder='/Users/hhunt1/Documents/Write ups/Maths/';
end
v = VideoWriter(strcat(savefolder,'Nucleus-cut',num2str(D),'-angle-',num2str(cutangle),'-wz.mp4'));
open(v)
figure
for t=2:size(cutca,2)
    plot([-h(end:-1:2) h],cutca{t},'LineWidth',2)
    xlabel('space (\mum)')
    ylabel('[Ca^{2+}] (\muM)')
    title(strcat('t=',num2str((t-1)*tint+tmin),'s'))
    axis([-max(h) max(h) min(min([cutca{2:end}])) max(max([cutca{2:end}]))])
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
end
function makemovie(tmin,tint,D,tsol)
if exist('/home/hhunt1/AnaNuc/1Dsims/data/','dir')
    savefolder='/home/hhunt1/AnaNuc/1Dsims/data/';
else
    savefolder='/Users/hhunt1/Documents/Write ups/Maths/';
end
v = VideoWriter(strcat(savefolder,'Nucleus-cyl',num2str(D),'-wz_cstq.avi'));
open(v)
figure
for t=1:size(tsol,2)
    surf(tsol{t},...
        'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(0,90)
    colorbar
    caxis([0.1 max(max([tsol{:}]))])
    title(strcat('t=',num2str((t-1)*tint+tmin),'s'))
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
end

function plotq(fits,crtpts)
xpts=cell(1,4);
figure
plot([1,crtpts(1)],fits{4}(crtpts(5))*ones(1,2))
hold on
for ii=1:4
    xpts{ii}=crtpts(ii):crtpts(ii+1);
    plot(xpts{ii},fits{ii}(xpts{ii}))
end
end
