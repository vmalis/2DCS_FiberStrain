%=========================================================================
%   Step2 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step1_2DSTR)
%
%           1)  mri_data.mat
%           2)  force_data.mat 
%           4)  strain rate eigen values from ev_calculated_data.mat
%           6)  fiber data
%           7)* optional progress file allowing you to strat from the 
%               folder where you've ended up priveously
%
% OUTput:   1) result.mat
%           2) roi.avi
%           3) image.png
%           4) Mask.mat
%_____________________________________________________
% required subroutines:
%       1)imoverlay
%       2)track2dv3
%       3)rdir
%_____________________________________________________
% 
% v.2.0 tailored for ULLS data analysis
% (improved from v.1.0 used for Age-related Differences in Strain Rate Tensor paper (MRM 2014))
% 
% written by Vadim Malis
% 01/15 at UCSD RIL
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%


% Number of rois per subject !PREDEFINED!
% below you can chanage the size of the rois go to line 
ROIs = 3;

switcher=true;

% Select major data folder
window_info=sprintf('Choose folder of the study');
RootPath = uigetdir('~/Desktop',window_info);
cd(RootPath)


% Uses progress file if not run for the first time
if exist(fullfile(cd, 'roi_progress.mat'), 'file') == 2
   load('roi_progress.mat')
   load(fiberfile)
else
    iii=1;
    
    %load fiber .mat file (output of step1)
    [fiberfile,fiberpath]=uigetfile('*.mat','select fiber data file');
    fiberfile=fullfile(fiberpath,fiberfile);
    file_list=rdir('**/mri_data.mat');
end     

% counter to check if repeat of roi placement or not
jj=0;

while switcher


PathName = file_list(iii).name;
PathName = PathName(1:end-13);
cd(RootPath)
cd(PathName)


if jj==0    % do only if not a repeat

    load('mri_data.mat')
    load('force_data.mat')
    load('ev_calculated_data.mat','EV_Negative')
    load('ev_calculated_data.mat','EV_Positive')
    load('ev_calculated_data.mat','EV_Sum')
    load('ev_calculated_data.mat','SR_Angle')
    load('ev_calculated_data.mat','v_ap_sm')
    load('ev_calculated_data.mat','v_rl_sm')
    load('ev_calculated_data.mat','v_si_sm')
    
    %contrast adjustment for magnitude images
    image=im_m;
    imgMin = double(min(image(:)));
    imgMax = double(max(image(:)));
    image  = (image - imgMin) / (imgMax - imgMin);
    
    % Loading data from ev_calculated_data
    %-------------------------------------
    angle = SR_Angle;
    %angles for preview
    ANGL = zeros(ROIs,numphases);

    nev = EV_Negative;
    %nev for preview
    NEV = zeros(ROIs,numphases);
    
    pev = EV_Positive;
    %PEV
    PEV = zeros(ROIs,numphases);

    sev = EV_Sum;
    %SEV
    SEV = zeros(ROIs,numphases);
    %-------------------------------------
    
    

%% Fibers check-adjust
check=true;


load(fiberfile)

%choose apropriate fiber data and perform tracking
sidx=find(strcmpi(PatientID,{fiber.id}));
slidx=find(abs([fiber(sidx).slice_location]-SliceLocation)<0.1);
if isempty(slidx)
fiber(size(fiber,2)+1).id=PatientID;
fiber(size(fiber,2)).slice_location=SliceLocation;
fiber(size(fiber,2)).x=128*ones(7,2);
fiber(size(fiber,2)).y=128*ones(7,2);
sidx=cat(2,sidx,size(fiber,2));
slidx=size(sidx,2);
end
fiber_x_0=round(fiber(sidx(slidx)).x/2);
fiber_y_0=round(fiber(sidx(slidx)).y/2);

while check

    [fiber_X,fiber_Y,~,~,~,~]=track2dv3(fiber_x_0(:),fiber_y_0(:),v_rl*-1,v_si,v_ap,dt,resolution);

    fiber_X=reshape(fiber_X,[size(fiber_X,1)/2,2,numphases]);
    fiber_Y=reshape(fiber_Y,[size(fiber_Y,1)/2,2,numphases]);       
    
    figure
    imshow(image(:,:,1),'InitialMagnification', 250);

    for j=1:numphases
        
        imshow(image(:,:,j),'Border','tight','InitialMagnification', 250);
        hold on
        for k=1:size(fiber_X,1)
            plot(fiber_X(k,:,j),fiber_Y(k,:,j),'-r');
            f_label=sprintf('%d',k);
            text(fiber_X(k,1,j)-15,fiber_Y(k,1,j),['\rightarrow' f_label],'Color','r','FontSize',16)
        end
        hold off
        pause(0.1)
        title('Fibers')
    end

    %questdlg
    choice = questdlg('Fibers endpoints are ok?','?', 'No','Yes','Yes');

    switch choice
        case 'No'
           
            close
            
            figure
            imshow(image(:,:,1),'Border','tight','InitialMagnification', 250)
            str1 = '\rightarrow ';
            hold on
                for k=1:size(fiber_X,1)
                    plot(fiber_X(k,:,1),fiber_Y(k,:,1),'-r');
                    str2 = num2str(k);
                    text(fiber_X(k,1,1)-15,fiber_Y(k,1,1),[str2 ' ' str1],'Color','r','FontSize',16);
                end
            
            f_n = inputdlg('Fibers you want to reposes?', '?');
            f_n=cell2mat(textscan(f_n{1},'%d','Delimiter',','));
            
            for q=1:size(f_n,1)
                   
            new_fiber_line=imline(gca);
            new_f=getPosition(new_fiber_line);
            
            fiber_x_0(f_n(q),1,1)=new_f(1,1);
            fiber_x_0(f_n(q),2,1)=new_f(2,1);
            fiber_y_0(f_n(q),1,1)=new_f(1,2);
            fiber_y_0(f_n(q),2,1)=new_f(2,2);
            
            end
           
            
        case 'Yes'
            
               %fiber angle
                Fiber_ANGL = zeros(ROIs,numphases);

                for kk=1:numphases
    
                    Fiber_ANGL(1,kk)=mean(fiberangle(fiber_X(1:2,1,kk),fiber_Y(1:2,1,kk),fiber_X(1:2,2,kk),fiber_Y(1:2,2,kk)));
                    Fiber_ANGL(2,kk)=mean(fiberangle(fiber_X(3:5,1,kk),fiber_Y(3:5,1,kk),fiber_X(3:5,2,kk),fiber_Y(3:5,2,kk)));
                    Fiber_ANGL(3,kk)=mean(fiberangle(fiber_X(6:7,1,kk),fiber_Y(6:7,1,kk),fiber_X(6:7,2,kk),fiber_Y(6:7,2,kk)));

                end
                
                check=false;

   end

close

end

fiber(sidx(slidx)).x=squeeze(fiber_X(:,:,1))*2;
fiber(sidx(slidx)).y=squeeze(fiber_Y(:,:,1))*2;

save(fiberfile,'fiber');


else

end

%---------fiber adjustment part ends here----------------------------------


%% Selecting ROIs

%roi confirmation variable
check=true;
while check

    
%%% Masking GUI 

%image stack
im_seq      =   zeros(r,c,numphases,3);
    
%converting to RGB
im_seq(:,:,:,1) = image;
im_seq(:,:,:,2) = image;
im_seq(:,:,:,3) = image;

I=squeeze(im_seq(:,:,1,:));
   


%-------------------------GUI----------------------------------------------
figure

imshow(I,'InitialMagnification', 200,'Border','tight');
hold on
for k=1:size(fiber_X,1)
       plot(fiber_X(k,:,1),fiber_Y(k,:,1),'-r');
end
hold off

title('1st magnitude image')


%%%------------ROI1--------------------------

%YOU CAN CHANGE ROI SIZE HERE----------------
roi_x=7;
roi_y=7;
%--------------------------------------------
h = imrect(gca,[100,100,roi_x,roi_y]);
setResizable(h,0)
wait(h);
mask   = createMask(h);
border = uint8(bwperim(mask,8));
[row1,col1] = find(mask);
I=imoverlay(I,border,[1 1 0]);

cla(gca)
imshow(I,'InitialMagnification', 200,'Border','tight');
hold on
for k=1:size(fiber_X,1)
       plot(fiber_X(k,:,1),fiber_Y(k,:,1),'-r');
end
hold off

title('1st magnitude image')


%%%------------ROI2--------------------------

%YOU CAN CHANGE ROI SIZE HERE----------------
roi_x=7;
roi_y=7;
%--------------------------------------------
h = imrect(gca,[col1(1),row1(1),roi_x,roi_y]);
setResizable(h,0)
wait(h);
mask   = createMask(h);
border = uint8(bwperim(mask,8));
[row2,col2] = find(mask);
I=imoverlay(I,border,[1 1 0]);

cla(gca)
imshow(I,'InitialMagnification', 200,'Border','tight');
hold on
for k=1:size(fiber_X,1)
       plot(fiber_X(k,:,1),fiber_Y(k,:,1),'-r');
end
hold off

title('1st magnitude image')

%%%------------ROI3--------------------------

%YOU CAN CHANGE ROI SIZE HERE----------------
roi_x=5;
roi_y=10;
%--------------------------------------------
h = imrect(gca,[col2(1),row2(1),roi_x,roi_y]);
setResizable(h,0)
wait(h);
mask   = createMask(h);
border = uint8(bwperim(mask,8));
[row3,col3] = find(mask);
I=imoverlay(I,border,[1 1 0]);

pause(0.1)
close


%---------GUI END----------------------------------------------------------






%tracking regions of interest
msg=msgbox('Please wait while tracking is in progress: ROI 1');
[xs1,ys1,~,~,~,~] = track2dv3(col1,row1,v_rl*-1,v_si,v_ap,dt,resolution);
close(msg)
msg=msgbox('Please wait while tracking is in progress: ROI 2');
[xs2,ys2,~,~,~,~] = track2dv3(col2,row2,v_rl*-1,v_si,v_ap,dt,resolution);
close(msg)
msg=msgbox('Please wait while tracking is in progress: ROI 3');
[xs3,ys3,~,~,~,~] = track2dv3(col3,row3,v_rl*-1,v_si,v_ap,dt,resolution);
close(msg)

%Creating MASK

MASK=zeros(r,c,numphases,ROIs);

for i=1:numphases
 
    
    for jj=1:size(row1,1)
        
    MASK(int16(ys1(jj,i)),int16(xs1(jj,i)),i,1)=1;
    im_seq(int16(ys1(jj,i)),int16(xs1(jj,i)),i,1)=1;
    im_seq(int16(ys1(jj,i)),int16(xs1(jj,i)),i,2)=1;
    im_seq(int16(ys1(jj,i)),int16(xs1(jj,i)),i,3)=0;

    end

    for jj=1:size(row2,1)
        
    MASK(int16(ys2(jj,i)),int16(xs2(jj,i)),i,2)=1;
    im_seq(int16(ys2(jj,i)),int16(xs2(jj,i)),i,1)=1;
    im_seq(int16(ys2(jj,i)),int16(xs2(jj,i)),i,2)=1;
    im_seq(int16(ys2(jj,i)),int16(xs2(jj,i)),i,3)=0;

    end
    
    
    for jj=1:size(row3,1)
        
    MASK(int16(ys3(jj,i)),int16(xs3(jj,i)),i,3)=1;
    im_seq(int16(ys3(jj,i)),int16(xs3(jj,i)),i,1)=1;
    im_seq(int16(ys3(jj,i)),int16(xs3(jj,i)),i,2)=1;
    im_seq(int16(ys3(jj,i)),int16(xs3(jj,i)),i,3)=0;

    end

end


%angle, nev, pev

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


sv  = zeros(ROIs,r,c,numphases);
for i=1:ROIs            %rois
    for j=1:numphases   %frames
    sv (i,:,:,j) =   sev(:,:,j).*MASK(:,:,j,i);
    end
end

sv (sv==0)  = NaN;
SEV(:,:) = squeeze(nanmean(nanmean(sv,2),3));
close


%% video

vid = VideoWriter('roi.avi');
vid.Quality=100;
open(vid);

FigHandle = figure('Position', [200, 100, 1000, 625]);
    for j=1:numphases
        
        subplot(3,3,[1 2 4 5 7 8]), imshow(squeeze(im_seq(:,:,j,:)),'Border','tight','InitialMagnification', 200);
        hold on
        for k=1:size(fiber_X,1)
        subplot(3,3,[1 2 4 5 7 8]), plot(fiber_X(k,:,j),fiber_Y(k,:,j),'-r');
        end
        title('Traced ROIs')
        vid_frame = getframe;

        
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
       
        writeVideo(vid,vid_frame);

    end

close(vid);

%------------------ questdlg----------------------------------------------
    choice = questdlg('ROIs placement is ok?','?', 'No','Yes','Yes');

    switch choice
        case 'No'
            check=true;
            jj=1;
        case 'Yes'
            check=false;
        name_image = sprintf('image.png'); 
        imwrite(I,name_image,'png');
        iii=iii+1;
        jj=0;
    end

close
end

%% questdlg quality grading will be used for averaging accros the sets
choice = questdlg('Image quality?','?', 'Poor','Good','Outstanding!',...
    'Outstanding!');

    switch choice
        case 'Poor'
            quality='poor';
        case 'Good'
            quality='good';  
        case 'Outstanding!'
            quality='outstanding';
    end

%% questdlg
choice = questdlg('More masks?','?', 'No','Yes','Yes');

    switch choice
        case 'No'
            switcher=false;
        case 'Yes'
            switcher=true;
    end
    
if iii>size(file_list,1)    
 switcher=false;
 h = msgbox('All the data are processed', 'Info','warn');
end
    
%% peak force permutation

ANGL=circshift(ANGL,numphases/2-frame,2);
NEV=circshift(NEV,numphases/2-frame,2);
PEV=circshift(PEV,numphases/2-frame,2);
SEV=circshift(SEV,numphases/2-frame,2);
Fiber_ANGL=circshift(Fiber_ANGL,numphases/2-frame,2);
 
%% saving
save('Mask.mat','MASK');
save('results.mat','PatientID');
save('results.mat','Series_name','-append');
save('results.mat','SliceLocation','-append');
save('results.mat','quality','-append');
save('results.mat','ANGL','-append');
save('results.mat','Fiber_ANGL','-append');
save('results.mat','NEV','-append');
save('results.mat','PEV','-append');  
save('results.mat','SEV','-append');  
  
    
cd ..
cd ..
cd ..
cd ..
cd ..

if exist(fullfile(cd, 'roi_progress.mat'), 'file') == 2
   save ('roi_progress.mat','file_list','-append')
   save ('roi_progress.mat','iii','-append')
else

save ('roi_progress.mat','file_list')
save ('roi_progress.mat','iii','-append') 
save ('roi_progress.mat','fiberfile','-append') 
end

end

if iii==size(file_list,1)
msg=msgbox('Analysis Completed!');
else
end

clear all
