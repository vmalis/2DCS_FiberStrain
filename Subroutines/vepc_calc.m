function [v_x,v_y,v_z]=vepc_calc(Data_phase,magnitude)

% velocity calculation from phase data
%   
% Input:   
%      phase data
%      magnitude data
%
% Output:  
%      arrays of reconstructed velocities
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis


%% create phase arrays (take only usefull voxels)
phase0=squeeze(Data_phase(:,:,:,1));
phasex=squeeze(Data_phase(:,:,:,2));
phasey=squeeze(Data_phase(:,:,:,3));
phasez=squeeze(Data_phase(:,:,:,4));

phase0(magnitude<.1)=0;
phasex(magnitude<.1)=0;
phasey(magnitude<.1)=0;
phasez(magnitude<.1)=0;


%% create velocity data

venc=10;

v_x=phasex-phase0;
v_y=phasey-phase0;
v_z=phasez-phase0;


%% unwrap
ix_max=find(v_x>pi);
ix_min=find(v_x<-pi);
iy_max=find(v_y>pi);
iy_min=find(v_y<-pi);
iz_max=find(v_z>pi);
iz_min=find(v_z<-pi);

v_x(ix_max)=v_x(ix_max)-2*pi;
v_x(ix_min)=v_x(ix_min)+2*pi;
v_y(iy_max)=v_y(iy_max)-2*pi;
v_y(iy_min)=v_y(iy_min)+2*pi;
v_z(iz_max)=v_z(iz_max)-2*pi;
v_z(iz_min)=v_z(iz_min)+2*pi;

%% convert to velocity
v_x=v_x/pi*venc;
v_y=v_y/pi*venc;
v_z=v_z/pi*venc;

% v_x=permute(v_x,[2,3,1]);
% v_y=permute(v_y,[2,3,1]);
% v_z=permute(v_z,[2,3,1]);


% correct for eddy current
error_x = im_mean(v_x);
error_y = im_mean(v_y);
error_z = im_mean(v_z);

v_x = im_subtract(v_x, error_x);
v_y = im_subtract(v_y, error_y);
v_z = im_subtract(v_z, error_z);



% 
% 
% figure; montage(mat2gray(abs(v_x)))
% figure; montage(mat2gray(abs(v_y)))
% figure; montage(mat2gray(abs(v_z)))