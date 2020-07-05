% H Hunt 2019
% Numerically solve buffered diffusion (nonlinear)
% Units used are ms, uM, um
% To do - check convergence with q_spline.
%           double check bd
function [sol,r,t]=nonlinearDiff_longer(D_uin,params_in,length_param,data_fits,data_crtpts,tint,tmax)
if exist('/Users/hhunt1/Documents/Nucleus/cutNuc1DData/ana_code','dir')
addpath('/Users/hhunt1/Documents/Nucleus/cutNuc1DData/ana_code')
end
global fits; global crtpts;global params; global D_u;global inclength;
fits=data_fits;
crtpts=data_crtpts;
inclength=length_param;
% params = [D_u,b_t,K,NPC,a];
if isempty(params_in)
    params_in=[1,0,1,1,2];
end
params=params_in;
D_u=D_uin;
m=1;
% a=4.7; % Av nucleus diameter measured by G Bass
% tmax=800;
% tsteps=201;%tmax/1+1;
r=linspace(0,params(5),2e2+1);
t=linspace(0,tmax,round(tmax/tint)+1);
sol=pdepe(m,@bufferDiff,@diffIC,@diffBC,r,t);
% sol = pdepe(m,@heatcyl,@heatic,@heatbc,x,t);
u=sol(:,:,1);
if 0
    figure
    mirror=1;
    if ~mirror
        surf(r(1:size(u,2)),t(1:size(u,1)),u,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');colorbar
    else
        surf([-r((size(u,2)):-1:2) r(1:size(u,2))],t(1:size(u,1)),[u(:,end:-1:2) u(:,1:end)],'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');colorbar
    end
    xlabel('x')
    ylabel('t')
    zlabel('c(x,t)')
    view([0,90])

%     figure
%     plot(t,q_spline(t,fits,crtpts))
%     figure
%     plot(t,sol(:,1))
%     xlabel('Time')
%     ylabel('calcium')
%     title('Concentration change at center of disc')
%     makemovie1D(0,t(2),D_u,u,r,0)
    plotq(fits,crtpts,u,t,1,D_u)
end
end
function [c,f,s]=bufferDiff(r,t,u,dudr)
global params; global D_u
b_t=params(2);%1;
K=params(3);
D_b=params(1);
theta=b_t*K./(K+u).^2;
c=(1+theta);
f=D_u*dudr+D_b*theta*dudr;
s=0;%-(t-200).^2/4e4+1;%q_spline(t);
end
function c0=diffIC(r)
global fits; global crtpts;
c0=q_spline(0,fits,crtpts);%0;%.1;
end
function [pl,ql,pr,qr]=diffBC(xl,ul,xr,ur,t)
global fits; global crtpts; global params; global inclength;
pl=0*ones(size(t));
ql=0*ones(size(t));
qr=0*ones(size(t));
ptime=fits{1}.breaks(7);
tprepeak=t<ptime;
tpostpeak=t>=ptime;
pr=t;
pr(tprepeak)=ur-params(4)*q_spline(t(tprepeak),fits,crtpts);%ur+(t-200).^2/4e4-1;
pr(tpostpeak)=ur-params(4)*q_spline(ptime+(t(tpostpeak)-ptime)*inclength,fits,crtpts);
end
function makemovie1D(tmin,tint,D,cutca,h,cutangle)
if exist('/home/hhunt1/AnaNuc/1Dsims/data/','dir')
    savefolder='/home/hhunt1/AnaNuc/1Dsims/data/';
else
%     savefolder='/Users/hhunt1/Documents/Write ups/Maths/';
    savefolder='/Users/hhunt1/Documents/huia/';
end
v = VideoWriter(strcat(savefolder,'Nucleus-num',num2str(D),'.avi'));
open(v)
figure
for t=2:size(cutca,1)
    xvals=h;%[-h(end:-1:2) h]
    plot(xvals,cutca(t,:),'LineWidth',2)
    xlabel('space (\mum)')
    ylabel('[Ca^{2+}] (\muM)')
    title(strcat('t=',num2str((t-1)*tint+tmin),'ms'))
    miny=min(cutca(:));
    axis([min(xvals) max(xvals) (miny-0.1*abs(miny)) 1.1*max(cutca(:))])
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)
end

function plotq(fits,crtpts,cutca,t,saveBool,D)
xpts=cell(1,3);
figure
qcolor=[0.8,0 0.5];
hold on
centreAv=mean(cutca(:,round((0.9*end/2):(1.1*end/2))),2);
h1=plot(t(1:size(centreAv,1)),centreAv,'Color',[qcolor(2:3),qcolor(1)],'LineWidth',2);
h2=plot(t,q_spline(t,fits,crtpts));
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
    if exist(savefolder)
    saveas(gcf,strcat(savefolder,'catrans_cn_D-',num2str(D),'-num.png'));
    end
end
end
