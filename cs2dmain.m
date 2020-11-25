%% -----------------------------------------------------------------------
%  UC San Diego / November 2020 / Vadim Malis
%
%   ISMRM 2021 Strain 2d analysis
%
%-------------------------------------------------------------------------

clc
clear all

disp("Starting reconstruction")
disp("=======================================")

%% read pfile
disp("Step 1: reading pfile to sturcture")
disp("========================================")
raw=pfile2struct;


%% recon 

Recon = reconCS(raw(1),'full',1);

for i=2:size(raw,2)
    
    ReconCS(i) = recon(raw(i), 'cs4vps2', i);
    
end

save cs_processed.mat Recon


%% read force
% USE: forceDataRead.m

force_data = forceDataRead(cd);

%MVC

%per series


%% read dicom image
%     USE:    dicom2struct.m

    % find corresponding slice based on location information
    
    % register highres image to dynamic based on location information
    % USE:    mri_register.m
   
    
    % plot fibers
    % get the fiber coordinates
    % plot on highres image, register highres to 1st frame of magnitude
    % dynamic and get fiber coordinates
    
    
    % track fibers
    % USE: track2dv4
    
    %rest look into code        /Subroutines2D/Step2_ROIs.mat
    % lines                 120 - 220

    
%% Strain analysis
% note in the codes it's done for 3d need 2d (no Z)

% USE:              Step03_Strain_calcs_v2
% lines             159-206

% strain calculations
% USE:             strain_tensor_calcs

