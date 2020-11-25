%=========================================================================
%   Step2 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strane Rate Toolkit
%
%               This is a kight version for "manual" preview-analysis
%
%=========================================================================
%
% The script performs 2d Strain Rate analysis per age group and also saves
% per subject information. Due to huge ammount of data to handle this
% script allows one-by-one mask placing so you wount get bored =)
% Once you are boared chose stop and continue next time, the script will
% check if the subject was already processed before
%
%
% INput: (output of Step1_2DSTR)
%
%           1) sorted velocities stack
%           2) filtered velocities stack
%
% OUTput:   1) avaraged per roi NEV, PEV, SR Angle
%
%_____________________________________________________
% required subroutines:
%       1)imoverlay
%       2)track2dv3
%_____________________________________________________
% 
% written by Vadim Malis
% 12/14 at UCSD RIL
%==========================================================================



% Number of rois per subject !PREDEFINED! Wount be changed
% below you can
ROIs = 3;


window_info=sprintf('Choose study folder');
PathName = uigetdir('~/Desktop',window_info);
cd(PathName)

load('mri_data.mat')

%angles
fname=sprintf('angle.dat');
fid = fopen(fname,'r','l');
angle = fread(fid,'float');
fclose(fid);
angle = reshape(angle, [r c numphases]);
angle = permute(angle, [2,1,3]);
ANGL = zeros(ROIs,numphases);


%NEV
fname=sprintf('negative_eig.dat');
fid = fopen(fname,'r','l');
nev = fread(fid,'float');
fclose(fid);
nev = reshape(nev, [r c numphases]);
nev = permute(nev, [2,1,3]);
NEV = zeros(ROIs,numphases);

%PEV
fname=sprintf('positive_eig.dat');
fid = fopen(fname,'r','l');
pev = fread(fid,'float');
fclose(fid);
pev = reshape(pev, [r c numphases]);
pev = permute(pev, [2,1,3]);
PEV = zeros(ROIs,numphases);






%% Selecting ROIs

%roi confirmation variable
check=true;
while check

    
%% Masking GUI

%contrast adjustment
image=im_m;
imgMin = double(min(image(:)));
imgMax = double(max(image(:)));
image  = (image - imgMin) / (imgMax - imgMin);  

%----------------------------------------
%image stack
im_seq      =   zeros(r,c,numphases,3);
    
%converting to RGB
im_seq(:,:,:,1) = image;
im_seq(:,:,:,2) = image;
im_seq(:,:,:,3) = image;

I=squeeze(im_seq(:,:,1,:));
    
figure
imshow(I,'InitialMagnification', 200);
title('1st magnitude image')

%%%   ROI1

%YOU CAN CHANGE ROI SIZE HERE----------------
h = imrect(gca,[128 128 7 7]);        
%-------------------------------------------- 

setResizable(h,0)
wait(h);
mask   = createMask(h);
border = uint8(bwperim(mask,8));
[row1,col1] = find(mask);
I=imoverlay(I,border,[1 1 0]);

close

%-------------------------------------------- 
%-------------------------------------------- 

figure
imshow(I,'InitialMagnification', 200);
title('1st magnitude image')

%%%   ROI2

%YOU CAN CHANGE ROI SIZE HERE----------------
h = imrect(gca,[128 128 7 7]);        
%-------------------------------------------- 

setResizable(h,0)
wait(h);
mask   = createMask(h);
border = uint8(bwperim(mask,8));
[row2,col2] = find(mask);
I=imoverlay(I,border,[1 1 0]);

close

%-------------------------------------------- 
%-------------------------------------------- 

figure
imshow(I,'InitialMagnification', 200);
title('1st magnitude image')

%%%   ROI3

%YOU CAN CHANGE ROI SIZE HERE----------------
h = imrect(gca,[128 128 5 10]);        
%-------------------------------------------- 

setResizable(h,0)
wait(h);
mask   = createMask(h);
border = uint8(bwperim(mask,8));
[row3,col3] = find(mask);
I=imoverlay(I,border,[1 1 0]);

msg=msgbox('Please wait while tracking is in progress');



%tracking
[xs1,ys1,vr,vx,vy,vz] = track2dv3(col1,row1,v_rl*-1,v_si,v_ap,dt,resolution);
[xs2,ys2,vr,vx,vy,vz] = track2dv3(col2,row2,v_rl*-1,v_si,v_ap,dt,resolution);
[xs3,ys3,vr,vx,vy,vz] = track2dv3(col3,row3,v_rl*-1,v_si,v_ap,dt,resolution);



%Creating MASK

MASK=zeros(r,c,numphases,ROIs);

for i=1:numphases
    
MASK(int16(ys1(:,i)),int16(xs1(:,i)),i,1)=1;
MASK(int16(ys2(:,i)),int16(xs2(:,i)),i,2)=1;
MASK(int16(ys3(:,i)),int16(xs3(:,i)),i,3)=1;

im_seq(int16(ys1(:,i)),int16(xs1(:,i)),i,1)=1;
im_seq(int16(ys1(:,i)),int16(xs1(:,i)),i,2)=1;
im_seq(int16(ys1(:,i)),int16(xs1(:,i)),i,3)=0;

im_seq(int16(ys2(:,i)),int16(xs2(:,i)),i,1)=1;
im_seq(int16(ys2(:,i)),int16(xs2(:,i)),i,2)=1;
im_seq(int16(ys2(:,i)),int16(xs2(:,i)),i,3)=0;

im_seq(int16(ys3(:,i)),int16(xs3(:,i)),i,1)=1;
im_seq(int16(ys3(:,i)),int16(xs3(:,i)),i,2)=1;
im_seq(int16(ys3(:,i)),int16(xs3(:,i)),i,3)=0;


end


%angles, nev, pev for preview

ang  = zeros(ROIs,r,c,numphases);
for i=1:ROIs            %rois
    for j=1:numphases   %frames
    ang (i,:,:,j) =   angle(:,:,j).*MASK(:,:,j,i);
    end
end

ang (ang==0)  = NaN;
ANGL(:,:) = squeeze(nanmean(nanmean(ang,2),3));



nv  = zeros(ROIs,r,c,numphases);
for i=1:ROIs            %rois
    for j=1:numphases   %frames
    nv (i,:,:,j) =   nev(:,:,j).*MASK(:,:,j,i);
    end
end

nv (nv==0)  = NaN;
NEV(:,:) = squeeze(nanmean(nanmean(nv,2),3));

pv  = zeros(ROIs,r,c,numphases);
for i=1:ROIs            %rois
    for j=1:numphases   %frames
    pv (i,:,:,j) =   pev(:,:,j).*MASK(:,:,j,i);
    end
end

pv (pv==0)  = NaN;
PEV(:,:) = squeeze(nanmean(nanmean(pv,2),3));








close(msg)
close


%% video

vid = VideoWriter('roi.avi');
vid.Quality=100;
open(vid);

FigHandle = figure('Position', [100, 100, 1200, 625]);
    for j=1:numphases
        
        subplot(3,3,[1 2 4 5 7 8]), imshow(squeeze(im_seq(:,:,j,:)),'Border','tight','InitialMagnification', 200);
        title('Traced ROIs')
        frame = getframe;
        
        subplot(3,3,3),plot(1:j,squeeze(ANGL(1,1:j)),'-r'), xlim([0 23]),ylim([-90 90])
        xlabel('frame number')
        tit=sprintf('Angle ROI 1');
        ylabel(tit);
        hold on
        
        subplot(3,3,6),plot(1:j,squeeze(ANGL(2,1:j)),'-r'), xlim([0 23]),ylim([-90 90])
        xlabel('frame number')
        tit=sprintf('Angle ROI 2');
        ylabel(tit);
        hold on
        
        subplot(3,3,9),plot(1:j,squeeze(ANGL(3,1:j)),'-r'), xlim([0 23]),ylim([-90 90])
        xlabel('frame number')
        tit=sprintf('Angle ROI 3');
        ylabel(tit);
        hold on
       
        writeVideo(vid,frame);

    end

close(vid);

    %questdlg
    choice = questdlg('ROIs placement is ok?','?', 'No','Yes','Yes');

    switch choice
        case 'No'
            check=true;
        case 'Yes'
            check=false;
        name_image = sprintf('image.png'); 
        imwrite(I,name_image,'png');
        i=i+1;
    end

close

end

clearvars -except PatientID NEV PEV ANGL


RESULT=cat(1,ANGL,PEV,NEV);
RESULT=RESULT';