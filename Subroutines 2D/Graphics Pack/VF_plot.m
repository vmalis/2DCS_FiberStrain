%=========================================================================
%   Additional script for Strain Rate Analysis Toolbox
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step1_2DSTR)
%
%           1)  ev_calculated.mat (subject specific)
%           
% OUTput:   
%
%           1) set of .png images 
%              with SR vector field overlayed over magnitude image      
%_____________________________________________________
% required subroutines:
%       1)multiwaitbar
%       2)freezeColors
%       3)export_fig package
%       4)quiver_color
%_____________________________________________________
%
% written by Vadim Malis
% 02/15 at UCSD RIL
%
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%



[FileName,PathName] = uigetfile('*.mat','choose a mat file','~/Desktop');

cd(PathName)
load (FileName);

start_frame=60;
end_frame=numphases;

multiWaitbar('Generating plots...', 0, 'Color', 'g');

Lambda1=permute(Y_1,[2,1,3]);Lambda1(~Lambda1)=NaN;
Lambda2=permute(Y_2,[2,1,3]);Lambda2(~Lambda2)=NaN;

NEV_u=Lambda1.*squeeze(VectorF_red(:,:,:,1));
NEV_v=Lambda1.*squeeze(VectorF_red(:,:,:,2));
PEV_u=Lambda2.*squeeze(VectorF_blue(:,:,:,1));
PEV_v=Lambda2.*squeeze(VectorF_blue(:,:,:,2));


[~,z_nev]=cart2pol(NEV_u,NEV_v);
[~,z_pev]=cart2pol(PEV_u,PEV_v);

Z_nev_max=max(z_nev(:));
Z_pev_max=max(z_pev(:));
Z_nev_min=min(z_nev(:));
Z_pev_min=min(z_pev(:));

Z_max=max(Z_nev_max,Z_pev_max);
Z_min=min(Z_nev_min,Z_pev_min);

% clear z_nev z_pev


[X,Y]=meshgrid(1:size(im_m,1),1:size(im_m,1));

units_nev = texlabel('-10^3 sec^(-1)');
units_pev = texlabel('10^3 sec^(-1)');


mkdir('A_SRmaps')
cd('A_SRmaps')

%%---------------------------NEV-------------------------------------------

t=0;

h_fig=figure('Visible','off');

for i=start_frame:end_frame;
    t=t+1;
    imshow(mat2gray(im_m(:,:,i)));
    hold on
    freezeColors
%    %auto colorbar  
%    quiver_color(X,Y,NEV_u(:,:,i),NEV_v(:,:,i),Z_min:Z_max*.5,units_nev,1);
    quiver_color(X,Y,NEV_u(:,:,i),NEV_v(:,:,i),0.001:1700,units_nev,1);
    title('Strain Rate: Negative')
    export_fig(sprintf('NEV#%d.png', i),'-m10','-png');
    hold off
    unfreezeColors
    multiWaitbar('Generating plots...', t/2/(end_frame));
end

close(h_fig)


%%---------------------------PEV-------------------------------------------
h_fig=figure('Visible','off');

for i=start_frame:end_frame;
    t=t+1;
    imshow(mat2gray(im_m(:,:,i)));
    hold on
    freezeColors
%    %auto colorbar   
%    quiver_color(X,Y,PEV_u(:,:,i),PEV_v(:,:,i),Z_min:Z_max*.5,units_pev,1)   
    quiver_color(X,Y,PEV_u(:,:,i),PEV_v(:,:,i),0.001:1700,units_pev,1)
    title('Strain Rate: Positive')
    export_fig(sprintf('PEV#%d.png', i),'-m10','-png');
    hold off
    unfreezeColors
    multiWaitbar('Generating plots...', t/2/(end_frame));
end

close(h_fig)

multiWaitbar('Generating plots...', 'Close');  

