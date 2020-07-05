% H Hunt 2019
% Used in findDvals

function [cutca,h,sol,rfun,zfun,tfun]=sim_nucleus_cyl_wz_spline(D,l,a,cutangle,crtpts,fits,tint,tmax,plotBool)
interval=0.01;
r=0:interval:a;
zpos=0:interval:l;
% maxn=100;
tmin=0;
% tmax=1;
% tint=0.01;
% t=[tmin:2.6:962,964.6:tint:tmax];
t=tmin:tint:tmax;
% add padding
%t=[tmin*ones(1,crtpts(1)) t(1:(end-crtpts(1)))];
[rfun,zfun,tfun,sol]=generatecaconc(a,l,D,r,t,crtpts,fits);
[cutca,h]=cutnuc(cutangle,l,a,r,zpos,sol,t);
if plotBool
%     makemovie(tmin,tint,D,rsol,r,phi)
    u=cellfun(@(x)x(:,round(end/2)),sol,'UniformOutput',false);
    u=[u{:}];
    u=reshape(u,[size(sol,2),size(rfun,1)]);
    figure
    surf([-r((size(u,2)):-1:2) r(1:size(u,2))],t(1:size(u,1)),[u(:,end:-1:2) u(:,1:end)],'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');colorbar
    xlabel('x')
    ylabel('t')
    zlabel('c(x,t)')
    view([0,0])
    makemovie1D(tmin,tint,D,cutca,h,cutangle)
    plotq(fits,crtpts,cutca,t,1,D,cutangle)
%     plotqte(t,crtpts,fits)
end
end
function [rfun,zfun,just_tfun,sol]=generatecaconc(a,l,D,r,t,crtpts,fits)
if isempty(D)
    D=[];
end
if size(fits,1)>1
    warning('too many fits?')
    fits=fits(1,:);
    crtpts=crtpts(1,:);
end
maxt=1e3;
sol=cell(size(t));
z=0:0.01:l;
alphan=zeros(1,maxt+1);
for n=1:maxt
    alphan(n+1) = fzero(@(z) besselj(0, a*z), [alphan(n)+0.1/a, alphan(n)+5/a]);
end
rfun=besselj(0,r'*alphan(2:end))./besselj(1,alphan(2:end)*a)./(alphan(2:end));
zfun=sin(z'*(2*(0:(maxt-1))+1)*pi/l)./(2*(0:(maxt-1))+1);
sol=cellfun(@(x)zeros(size(rfun(:,1)*zfun(:,1)')),sol,'Un',0);
common_coeff=8/(pi*a);
just_tfun=cell(1,size(t,2));
for time=1:size(t,2)
    n=1:maxt;
    m=0:(maxt-1);
    intofe=int_e_fit(D,m,n,alphan,l,crtpts,fits,t(time));
    just_tfun{time}=-common_coeff*(intofe-0*exp(-D*(alphan(n+1)'.^2+((2*m+1)*pi/l).^2)*t(time)));
    sol{time}=rfun*just_tfun{time}*zfun';
    if sum(isnan(sol{time}(:)))>0
        warning('tsol isnan')
    end
    % Work out which fit applies
    fitNum=find(crtpts>t(time),1,'first')-1;
    if fitNum==0
        qadd=ppval(fits{1},fits{1}.breaks(1));
    elseif fitNum >1 &fitNum<3
        qadd=ppval(fits{1},(t(time)));
    elseif fitNum<=size(fits{1}.coefs,1)+1
        qadd=ppval(fits{1},(t(time)));
    else
        qadd=ppval(fits{1},(t(time)));
    end
    NE_release=0;%.1;%2*(qadd>0.1);
    sol{time}=sol{time}+qadd+NE_release;
    if sum(abs(sol{time}(:))>2*(qadd+NE_release))>0
        warning('something is not right with equations')
    end
    if sum(sol{time}(:)<0)>0
        warning(strcat('negative concentrations at timestep:',num2str(time)));
    end
%     tsol{time}=sum(cat(3,tfun{:}),3);
end
end

function sol=int_e_fit(D,m,n,alphan,l,crtpts,fits,time)
ecoeff=D*(alphan(n+1)'.^2+((2*m+1)*pi/l).^2);
% Work out how many segments to include
numSegs=find(crtpts>(time),1,'first');
if isempty(numSegs)
    numSegs=crtpts(end);
end
% v_0=fits{2}(crtpts(end));
if numSegs==1
    sol=0;%v_0*(1-exp(-ecoeff*time))./ecoeff;
else
% Calculate integral for each segment
try
segInts=zeros(size(n,2),size(m,2),numSegs);
catch
    warning('work out problem')
end
% segInts(:,:,1)=0;%v_0*(exp(ecoeff*(crtpts(1)-time))-exp(-ecoeff*time))./ecoeff;
% These fits are of the form p1x^2+p2x+p3
% We are taking the definite integral of the derivative*e^~
% int((2p1x+p2)e~)=...
for ii=2:(numSegs-1)
    try
%     segInts(:,:,ii)=2*fits{ii-1}.p1*((ecoeff*crtpts(ii)-1).*exp(ecoeff*(crtpts(ii)-time)) ...
%         -(ecoeff*crtpts(ii-1)-1).*exp(ecoeff*(crtpts(ii-1)-time)))./ecoeff.^2 ...
%     +fits{ii-1}.p2*(exp(ecoeff*(crtpts(ii)-time))-exp(ecoeff*(crtpts(ii-1)-time)))./ecoeff;
%     segInts(:,:,ii)=(exp(ecoeff*(crtpts(ii)-time))*(2*fits{ii-1}.p1*(ecoeff*crtpts(ii)-1)+fits{ii-1}.p2)...
%         -exp(ecoeff*(crtpts(ii-1)-time))*(2*fits{ii-1}.p1*(ecoeff*crtpts(ii-1)-1)+fits{ii-1}.p2))./ecoeff.^2;
    segInts(:,:,ii)=(exp(ecoeff*(crtpts(ii)-time)).*...
        (3*fits{1}.coefs(ii-1,1)*(ecoeff.^2*crtpts(ii)^2-2*ecoeff*crtpts(ii)+2)...
        +2*(fits{1}.coefs(ii-1,2)-3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1))*ecoeff.*(ecoeff*crtpts(ii)-1)...
        +ecoeff.^2*(fits{1}.coefs(ii-1,3)+3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1).^2-2*fits{1}.coefs(ii-1,2)*fits{1}.breaks(ii-1)))...
        -exp(ecoeff*(crtpts(ii-1)-time)).*...
        (3*fits{1}.coefs(ii-1,1)*(ecoeff.^2*crtpts(ii-1)^2-2*ecoeff*crtpts(ii-1)+2)...
        +2*(fits{1}.coefs(ii-1,2)-3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1))*ecoeff.*(ecoeff*crtpts(ii-1)-1)...
        +ecoeff.^2*(fits{1}.coefs(ii-1,3)+3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1).^2-2*fits{1}.coefs(ii-1,2)*fits{1}.breaks(ii-1))))./ecoeff.^3;
%         segInts(:,:,ii)=((3*fits{1}.coefs(ii-1,1)*(ecoeff.^2*crtpts(ii)^2-2*ecoeff*crtpts(ii)+2)...
%         +2*(fits{1}.coefs(ii-1,2)-3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1))*ecoeff.*(ecoeff*crtpts(ii)-1)...
%         +ecoeff.^2*(fits{1}.coefs(ii-1,3)+3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1).^2-2*fits{1}.coefs(ii-1,2)*fits{1}.breaks(ii-1)))...
%         -exp(ecoeff*(crtpts(ii-1)-crtpts(ii))).*...
%         (3*fits{1}.coefs(ii-1,1)*(ecoeff.^2*crtpts(ii-1)^2-2*ecoeff*crtpts(ii-1)+2)...
%         +2*(fits{1}.coefs(ii-1,2)-3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1))*ecoeff.*(ecoeff*crtpts(ii-1)-1)...
%         +ecoeff.^2*(fits{1}.coefs(ii-1,3)+3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1).^2-2*fits{1}.coefs(ii-1,2)*fits{1}.breaks(ii-1))))./ecoeff.^3;
    catch
        warning('weird fit?')
    end
end
ii=numSegs;
if 0%numSegs==(size(crtpts,2))
    segInts(:,:,ii)=fits{2}.a*fits{2}.b*...
        (exp(fits{2}.b*time)-...
        exp(ecoeff*(crtpts(ii-1)-time)+fits{2}.b*crtpts(ii-1)))./(ecoeff+fits{2}.b)...
        +fits{2}.c*fits{2}.d*...
        (exp(fits{2}.d*time)-...
        exp(ecoeff*(crtpts(ii-1)-time)+fits{2}.d*crtpts(ii-1)))./(ecoeff+fits{2}.d);
%     segInts(:,:,ii)=fits{ii-1}.a*fits{ii-1}.b*(exp(fits{ii-1}.b*crtpts(ii)+ecoeff*(crtpts(ii)-time))...
%         -exp(fits{ii-1}.b*crtpts(ii-1)+ecoeff*(crtpts(ii-1)-time)))./(ecoeff+fits{ii-1}.b)+...
%         fits{ii-1}.c*fits{ii-1}.d*(exp(fits{ii-1}.d*crtpts(ii)+ecoeff*(crtpts(ii)-time))...
%         -exp(fits{ii-1}.d*crtpts(ii-1)+ecoeff*(crtpts(ii-1)-time)))./(ecoeff+fits{ii-1}.d);
else
%     segInts(:,:,ii)=2*fits{ii-1}.p1*((ecoeff*crtpts(ii)-1) ...
%         -(ecoeff*crtpts(ii-1)-1).*exp(ecoeff*(crtpts(ii-1)-time)))./ecoeff.^2 ...
%     +fits{ii-1}.p2*(1-exp(ecoeff*(crtpts(ii-1)-time)))./ecoeff;
    try
    segInts(:,:,ii)=((3*fits{1}.coefs(ii-1,1)*(ecoeff.^2*time^2-2*ecoeff*time+2)...
        +2*(fits{1}.coefs(ii-1,2)-3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1))*ecoeff.*(ecoeff*time-1)...
        +ecoeff.^2*(fits{1}.coefs(ii-1,3)+3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1).^2-2*fits{1}.coefs(ii-1,2)*fits{1}.breaks(ii-1)))...
        -exp(ecoeff*(crtpts(ii-1)-time)).*...
        (3*fits{1}.coefs(ii-1,1)*(ecoeff.^2*crtpts(ii-1)^2-2*ecoeff*crtpts(ii-1)+2)...
        +2*(fits{1}.coefs(ii-1,2)-3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1))*ecoeff.*(ecoeff*crtpts(ii-1)-1)...
        +ecoeff.^2*(fits{1}.coefs(ii-1,3)+3*fits{1}.coefs(ii-1,1)*fits{1}.breaks(ii-1).^2-2*fits{1}.coefs(ii-1,2)*fits{1}.breaks(ii-1))))./ecoeff.^3;
    catch
        warning('weird fit?')
    end
end
if sum(isinf(segInts(:)))>0
    warning('something inf')
end
if time==171.5
sol=sum(segInts,3);
else
sol=sum(segInts,3);
end
end
end
function [sol,h]=cutnuc(cutangle,l,a,r,z,sol2d,t)
opponadj=tan(cutangle);
if opponadj<0
    cutangle=pi-mod(cutangle,pi);
    opponadj=tan(cutangle);
end
interval=0.01;
if opponadj==0
    hyplen=l;
    h=0:interval:hyplen;
    y=0*h;
    x=round(0:(l/(size(h,2)-1)/2):(l),2);
    idz=arrayfun(@(s)find(s==round(z,2),1,'first'),round(x,2),'UniformOutput',false);
    idz=[idz{:}];
    idr=arrayfun(@(s)find(s==round(r,2),1,'first'),round(y,2),'UniformOutput',false);
    idr=[idr{:}];
elseif 2*a/l<=opponadj
    warning('redo this section of if statement')
    hyplen=sqrt(a^2+(a/opponadj)^2);
    h=0:interval:hyplen;
%     y=round(0:(a/(size(h,2)-1)):a,2);
%     x=round(0:(l*opponadj/(4*(size(h,2)-1))):(l*opponadj/2),2);
    
    x=round(0:(a/(opponadj*(size(h,2)-1))):(2*a/(opponadj)),2); % l dir
    y=round(0:(a/(size(h,2)-1)):(a),2);
    idz=arrayfun(@(s)find(s==round(z,2),1,'first'),round(x,2),'UniformOutput',false);
    idz=[idz{:}];
    idr=arrayfun(@(s)find(s==round(r,2),1,'first'),round(y,2),'UniformOutput',false);
    idr=[idr{:}];
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
sol=cell(size(t));
if ~isempty(y)
    for time=1:size(t,2)
        try
%         sol{time}=tsol{time}(idx([end:-1:2 2 2:end]),idy)+integral(fun,0,t(time));
        sol{time}=sol2d{time}((idz-1)*size(sol2d{time},1)+(idr([end:-1:2 2 2:end])));%integral(fun,0,t(time));
%         sol{time}=rsol(idx([end:-1:2 2 2:end]))...
%             .*zsol(idy)*tsol{time}+integral(fun,0,t(time));
        catch
            warning('oh no!')
        end
    end
else
    for time=1:size(t,2)
        warning('add zcoeff at line 90,sim_nuc_cyl')
        sol{time}=sol2d{time}(2*ones(1,2*size(h,2)-1),ones(1,2*size(h,2)-1))+integral(fun,0,t(time));
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
v = VideoWriter(strcat(savefolder,'Nucleus-cut',num2str(D),'-angle-',num2str(cutangle),'-wz.avi'));
open(v)
figure
for t=2:size(cutca,2)
    plot([-h(end:-1:2) h],cutca{t},'LineWidth',2)
    xlabel('space (\mum)')
    ylabel('[Ca^{2+}] (\muM)')
    title(strcat('t=',num2str((t-1)*tint+tmin),'ms'))
    axis([-max(h) max(h) 0.9*min(min([cutca{:}])) 1.1*max(max([cutca{:}]))])
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
end
function makemovie(tmin,tint,D,caconc,r,phi)
if exist('/home/hhunt1/AnaNuc/1Dsims/data/','dir')
    savefolder='/home/hhunt1/AnaNuc/1Dsims/data/';
else
    savefolder='/Users/hhunt1/Documents/Write ups/Maths/';
end
v = VideoWriter(strcat(savefolder,'Nucleus-cyl',num2str(D),'-wz.avi'));
open(v)
figure
[R,PHI] = meshgrid(r,phi);
for t=1:size(caconc,2)
    Z=repmat(caconc{t},[size(phi,2),1]);
    surf(R.*cos(PHI), R.*sin(PHI), Z,...
        'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
    view(0,90)
    colorbar
    caxis([0.1 1.1])
    title(strcat('t=',num2str((t-1)*tint+tmin),'ms'))
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
end

function plotq(fits,crtpts,cutca,t,saveBool,D,cutangle)
xpts=cell(1,3);
figure
qcolor=[0.8,0 0.5];
hold on
centreAv=[cellfun(@(x)mean(x(800:1009)),cutca,'UniformOutput',false)];
h1=plot(t,[centreAv{:}],'Color',[qcolor(2:3),qcolor(1)],'LineWidth',2);
h2=plot([1,crtpts(1,1)],fits{1,2}(crtpts(1,end))*ones(1,2),'Color',qcolor,'LineWidth',2);
xpts{1}=crtpts(1,2):crtpts(1,end-1);
plot(xpts{1},ppval(fits{1,1},xpts{1}),'Color',qcolor,'LineWidth',2)
xpts{2}=crtpts(1,end-1):crtpts(1,end);
plot(xpts{2},fits{1,2}(xpts{2}),'Color',qcolor,'LineWidth',2)
xlabel('time (ms)','FontSize',18)
ylabel('[Ca^{2+}] (\muM)','FontSize',18)
legend([h2,h1],{'cytosol','nucleus'},'FontSize',18)
title(strcat('D=',num2str(round(D,2)),'\mum^2ms^{-1}'),'FontSize',24)
% xpts{3}=fits{1,1}.breaks;
% plot(xpts{3},ppval(fits{1,1},xpts{3}),'x')
if saveBool
    if exist('/home/hhunt1/AnaNuc/1Dsims/data/','dir')
        savefolder='/home/hhunt1/AnaNuc/1Dsims/data/';
    else
        savefolder='/Users/hhunt1/Documents/Write ups/Maths/';
    end
    saveas(gcf,strcat(savefolder,'catrans_cn_D-',num2str(D),'_angle-',num2str(cutangle),'-wz.png'));
end
end

function plotqte(t,crtpts,fits)
% ecoeff=D*(alphan(n+1)'.^2+((2*m+1)*pi/l).^2);
xpts=cell(1,3);
ypts=xpts;
timeStep=0.26;
% t=1:crtpts(end);
fit1der=fnder(fits{1},1);
figure
plot([t(1),crtpts(1,1)],fits{1,2}(crtpts(1,end))*zeros(1,2))
xpts{1}=fits{1}.breaks(1):0.01:fits{1}.breaks(end);
hold on

ypts{1}=ppval(fit1der,xpts{1}).*exp(-xpts{1});

xpts{2}=crtpts(1,end-1):0.01:crtpts(1,end);
ypts{2}=ppval(fnder(fits{2},1),xpts{2}).*exp(-xpts{2});

xpts{3}=fits{1}.breaks;
ypts{3}=ppval(fit1der,xpts{3}).*exp(-xpts{3});
for ii =1
    plot(xpts{ii},ypts{ii})
end
plot(xpts{3},ypts{3},'x')
xlim([t(1) t(end)])
end
