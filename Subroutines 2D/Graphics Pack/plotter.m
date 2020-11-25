%=========================================================================
%   Additional script for Strain Rate Analysis Toolbox
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step4_Analysis)
%
%           1)  Data_array.mat
%           
% OUTput:   
%           1) set of plots  
%_____________________________________________________
%
% written by Vadim Malis
% 12/15 at UCSD RIL
%
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%

%ecc
angl_lim=[30 100];
%iso
%angl_lim=[10 140];

%ecc
nev_lim=[-350 0];
%iso
%nev_lim=[-750 0];

%ecc
pev_lim=[0 350];
%iso
%pev_lim=[0 750];

%ecc
sev_lim=[-250 250];
%iso
%sev_lim=[-450 450];




nev=squeeze(nanmean(nanmean(NEV,1),3));
pev=squeeze(nanmean(nanmean(PEV,1),3));
sev=(-1)*squeeze(nanmean(nanmean(SEV,1),3));
ang=squeeze(nanmean(nanmean(FR_ANGL-SR_ANGL,1),3));

%std
nev_std=squeeze(nanstd(squeeze(nanmean(NEV,3)),1));
pev_std=squeeze(nanstd(squeeze(nanmean(PEV,3)),1));
sev_std=squeeze(nanstd(squeeze(nanmean(SEV,3)),1));
ang_std=squeeze(nanstd(squeeze(nanmean(FR_ANGL-SR_ANGL,3)),1));

%not including post data
nev(3,:,:)=[];
pev(3,:,:)=[];
sev(3,:,:)=[];
ang(3,:,:)=[];

nev_std(3,:,:)=[];
pev_std(3,:,:)=[];
sev_std(3,:,:)=[];
ang_std(3,:,:)=[];




% ploting each of three ROI on same plot...

hnv=figure;

intData1 = interp1(1:22,squeeze(nev(1,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(nev(1,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(nev(1,3,:)),1:0.1:22,'spline');

hold on

errorbar(squeeze(nev(1,1,:)),squeeze(nev_std(1,1,:))/3,'or');
errorbar(squeeze(nev(1,2,:)),squeeze(nev_std(1,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(nev(1,3,:)),squeeze(nev_std(1,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')

ylim(nev_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Negative Eigen Value');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','southeast');
set(gcf,'color','w');
filename=sprintf('NEV_pre.eps');
export_fig(filename,'-eps');
hold off



hnv=figure;

intData1 = interp1(1:22,squeeze(nev(2,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(nev(2,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(nev(2,3,:)),1:0.1:22,'spline');

hold on

errorbar(squeeze(nev(2,1,:)),squeeze(nev_std(2,1,:))/3,'or');
errorbar(squeeze(nev(2,2,:)),squeeze(nev_std(2,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(nev(2,3,:)),squeeze(nev_std(2,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')

ylim(nev_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Negative Eigen Value');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','southeast');
set(gcf,'color','w');
filename=sprintf('NEV_post.eps');
export_fig(filename,'-eps');
hold off


% ploting each of three ROI on same plot...

hpv=figure;

intData1 = interp1(1:22,squeeze(pev(1,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(pev(1,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(pev(1,3,:)),1:0.1:22,'spline');

hold on

errorbar(squeeze(pev(1,1,:)),squeeze(pev_std(1,1,:))/3,'or');
errorbar(squeeze(pev(1,2,:)),squeeze(pev_std(1,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(pev(1,3,:)),squeeze(pev_std(1,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')

ylim(pev_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Positive Eigen Value');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','northeast');
set(gcf,'color','w');
filename=sprintf('PEV_pre.eps');
export_fig(filename,'-eps');
hold off



hpv=figure;

intData1 = interp1(1:22,squeeze(pev(2,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(pev(2,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(pev(2,3,:)),1:0.1:22,'spline');

hold on

errorbar(squeeze(pev(2,1,:)),squeeze(pev_std(2,1,:))/3,'or');
errorbar(squeeze(pev(2,2,:)),squeeze(pev_std(2,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(pev(2,3,:)),squeeze(pev_std(2,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')

ylim(pev_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Positive Eigen Value');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','northeast');
set(gcf,'color','w');
filename=sprintf('PEV_post.eps');
export_fig(filename,'-eps');
hold off




% ploting each of three ROI on same plot...
 
hsv=figure;
 
intData1 = interp1(1:22,squeeze(sev(1,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(sev(1,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(sev(1,3,:)),1:0.1:22,'spline');
 
hold on
 
errorbar(squeeze(sev(1,1,:)),squeeze(sev_std(1,1,:))/3,'or');
errorbar(squeeze(sev(1,2,:)),squeeze(sev_std(1,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(sev(1,3,:)),squeeze(sev_std(1,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')
 
ylim(sev_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Sum Eigen Value');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','northeast');
set(gcf,'color','w');
filename=sprintf('sev_pre.eps');
export_fig(filename,'-eps');
hold off
 
 
 
hsv=figure;
 
intData1 = interp1(1:22,squeeze(sev(2,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(sev(2,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(sev(2,3,:)),1:0.1:22,'spline');
 
hold on
 
errorbar(squeeze(sev(2,1,:)),squeeze(sev_std(2,1,:))/3,'or');
errorbar(squeeze(sev(2,2,:)),squeeze(sev_std(2,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(sev(2,3,:)),squeeze(sev_std(2,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')
 
ylim(sev_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Sum Eigen Value');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','northeast');
set(gcf,'color','w');
filename=sprintf('sev_post.eps');
export_fig(filename,'-eps');
hold off
 


 
 
% ploting each of three ROI on same plot...
 
hag=figure;
 
intData1 = interp1(1:22,squeeze(ang(1,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(ang(1,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(ang(1,3,:)),1:0.1:22,'spline');
 
hold on
 
errorbar(squeeze(ang(1,1,:)),squeeze(ang_std(1,1,:))/3,'or');
errorbar(squeeze(ang(1,2,:)),squeeze(ang_std(1,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(ang(1,3,:)),squeeze(ang_std(1,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')
 
ylim(angl_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Strain Rate - Fiber Angle');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','northeast');
set(gcf,'color','w');
filename=sprintf('ANG_pre.eps');
export_fig(filename,'-eps');
hold off
 
 
 
hag=figure;
 
intData1 = interp1(1:22,squeeze(ang(2,1,:)),1:0.1:22,'spline');
intData2 = interp1(1:22,squeeze(ang(2,2,:)),1:0.1:22,'spline');
intData3 = interp1(1:22,squeeze(ang(2,3,:)),1:0.1:22,'spline');
 
hold on
 
errorbar(squeeze(ang(2,1,:)),squeeze(ang_std(2,1,:))/3,'or');
errorbar(squeeze(ang(2,2,:)),squeeze(ang_std(2,2,:))/3,'o','Color',[0,0.5,0]);
errorbar(squeeze(ang(2,3,:)),squeeze(ang_std(2,3,:))/3,'ob');
plot(1:0.1:22,intData1,'-r')
plot(1:0.1:22,intData2,'-','Color',[0,0.5,0])
plot(1:0.1:22,intData3,'-b')
 
ylim(angl_lim)
xlim([0 23])
xlabel('frame number')
tit=sprintf('Strain Rate - Fiber Angle');
ylabel(tit);
legend('Proximal','Middle','Distal','Location','northeast');
set(gcf,'color','w');
filename=sprintf('ANG_post.eps');
export_fig(filename,'-eps');
hold off
 


