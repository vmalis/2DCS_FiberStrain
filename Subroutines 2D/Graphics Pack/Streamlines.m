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
%           1) set of .png images 
%              with SR vector field overlayed over magnitude image      
%_____________________________________________________
% required subroutines:
%       1)multiwaitbar
%       2)freezeColors
%       3)export_fig package
%       4)StreamColor
%_____________________________________________________
%
% written by Vadim Malis
% 02/15 at UCSD RIL
%
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%

%frame number
N=2;


[FileName,PathName] = uigetfile('*.mat','choose a mat file');

cd(PathName)
load (FileName);
load ('mri_data.mat');

Y_1=permute(EV_Negative,[2,1,3]);
Y_2=permute(EV_Positive,[2,1,3]);
VectorF_red=NEVector;
VectorF_blue=PEVector;

% multiWaitbar('Generating plots...', 0, 'Color', 'g');

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

I=squeeze(im_m(:,:,N));


[X,Y]=meshgrid(1:size(im_m,1),1:size(im_m,1));

figure
imshow(mat2gray(I),'InitialMagnification',300);
%h=imellipse;
%h=impoly;

h=drawpolygon();
MASK=createMask(h);

MASK2=zeros(size(MASK));
% for i=1:size(MASK,2)
% j=find(MASK(i,:), 1 );
% MASK2(i,j:j+12)=1;
% end



index=find(MASK>0);
[Sy,Sx]=ind2sub([256,256],index);

Sy(1:2:end)=[];
Sx(1:2:end)=[];
Sy(1:2:end)=[];
Sx(1:2:end)=[];

imshow(mat2gray(I),'InitialMagnification',600)
colormap(gray)
freezeColors
colormap(jet)
hold on


Streamcolor(X,Y,squeeze(PEV_u(:,:,N)).*MASK,squeeze(PEV_v(:,:,N)).*MASK,Sx,Sy,squeeze(z_pev(:,:,N)).*MASK)

% hold on
%         for k=1:size(fiber_X,1)
%             plot(fiber_X(k,:,N),fiber_Y(k,:,N),'-k','LineWidth',.5);
%         end
% hold off

filename=sprintf('PEV%2d.png', N);


export_fig pev3 '-png' -m4

% close all
% 
% 
% imshow(mat2gray(I),'InitialMagnification',500)
% colormap(gray)
% freezeColors
% colormap(jet)
% hold on
% 
% Streamcolor(X,Y,squeeze(NEV_u(:,:,N)).*MASK,squeeze(NEV_v(:,:,N)).*MASK,Sx,Sy,squeeze(z_nev(:,:,N)).*MASK)
% filename=sprintf('NEV%2d.png', N);
% export_fig nev11 '-png' -m4
% 

close all