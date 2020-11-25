%=========================================================================
%   Additional script for Strain Rate Analysis Toolbox
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step2_ROIs)
%
%           1)  ev_calculated.mat (subject specific)
%           2)  overal_results.mat 
%               (to get info on max force frame and max eigenvalue along the fiber)
%
% OUTput:   1) vector field maps *.pdf or eps for spcified frame
%_____________________________________________________
% required subroutines:
%       1)rdir
%_____________________________________________________
% 
% v.2.0 tailored for ULLS data analysis
% (improved from v.1.0 used for Age-related Differences in Strain Rate Tensor paper (MRM 2014))
% 
% written by Vadim Malis
% 12/14 at UCSD RIL
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%


[FileName,PathName] = uigetfile('*.mat','choose a mat file','~/Desktop');

cd(PathName)
load (FileName);

start_frame=1;
end_frame=10;

%Select roi
figure
imshow(mat2gray(im_m(:,:,1)),'InitialMagnification', 200)
image_crop=imrect(gca);
image_crop_pos=getPosition(image_crop);

CROP_X=int16([image_crop_pos(1),image_crop_pos(1)+image_crop_pos(3)]);
CROP_Y=int16([image_crop_pos(2),image_crop_pos(2)+image_crop_pos(4)]);


close




multiWaitbar('Generating plots...', 0, 'Color', 'g');


% Open video files to write in
writerObj1 = VideoWriter('NEV_SL.avi');
% writerObj2 = VideoWriter('PEV_SL.avi');
open(writerObj1);
% open(writerObj2);


for i=start_frame:end_frame;
    for x=1:256
        for y=1:256
             Vector_quiv_u_red (x,y) = VectorF_red (x,y,i,1);
             Vector_quiv_v_red (x,y) = VectorF_red (x,y,i,2);
             Vector_quiv_u_blue(x,y) = VectorF_blue(x,y,i,1);
             Vector_quiv_v_blue(x,y) = VectorF_blue(x,y,i,2);
        end
    end
% h_fig1=figure('Visible','off');
% imshow(Y_1(:,:,i)','DisplayRange',[-1000 200],'InitialMagnification',250);
% colormap(winter)
% colorbar
% hold on
[X,Y] = meshgrid(1:256,1:256);
% quiver(X,Y,Vector_quiv_u_red,Vector_quiv_v_red,1,'r');
% axis tight
% axis([CROP_X(1) CROP_X(2) CROP_Y(1) CROP_Y(2)])
% hold off
% fullname=sprintf('NEV#%d.png',i);
% print(h_fig1, '-dpng', fullname);
% 
% h_fig2=figure('Visible','off');
% imshow(Y_2(:,:,i)','DisplayRange',[-200 1000],'InitialMagnification',250);
% colormap(autumn)
% colorbar
% hold on
% quiver(X,Y,Vector_quiv_u_blue,Vector_quiv_v_blue,1,'b');
% axis tight
% axis([CROP_X(1) CROP_X(2) CROP_Y(1) CROP_Y(2)])
% hold off
% fullname=sprintf('PEV#%d.png',i);
% print(h_fig2, '-dpng', fullname);
% 

[x,y] =meshgrid(CROP_X(1):4:CROP_X(2),CROP_Y(1):4:CROP_Y(2));

x=double(x);
y=double(y);

h_fig3=figure('Visible','off');
% imshow(mat2gray(im_m(:,:,i)),'InitialMagnification', 200)
hold on
hlines=streamcolor(X,Y,Vector_quiv_u_red,Vector_quiv_v_red,x,y,Y_2(:,:,i)');
colorbar
axis off
axis square
axis([CROP_X(1) CROP_X(2) CROP_Y(1) CROP_Y(2)])
% set(hlines,'LineWidth',1,'Color','r')
hold off
writeVideo(writerObj1,getframe(gca));

% h_fig4=figure('Visible','off');
% imshow(mat2gray(im_m(:,:,i)),'InitialMagnification', 200)
% hold on
% hlines=streamline(X,Y,Vector_quiv_u_blue,Vector_quiv_v_blue,x,y);
% set(hlines,'LineWidth',1,'Color','b')
% axis tight
% axis([CROP_X(1) CROP_X(2) CROP_Y(1) CROP_Y(2)])
% hold off
% writeVideo(writerObj2,getframe(gca));

multiWaitbar('creating plots and movies', i/(end_frame));

end

multiWaitbar('creating plots and movies', 'Close');  

% Close video files
close(writerObj1);
% close(writerObj2);
