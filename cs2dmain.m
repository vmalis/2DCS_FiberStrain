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


%% recon  data

% 1st set is full space
Recon = reconCS(raw(1),'full',1);

% sets 2-5 are CS
for i=2:size(raw,2)
    Recon_CS(i-1) = reconCS(raw(i), 'cs4vps2', i); 
end



%% anisotrpoic diffusion filter 
% from here for all calculation the smoothed velocity volumes must be used

num_iter=10; kappa=2; option=1; delta_t=1/7;
% calling the anisodiff2d code; options are set above; 
% less smoothing with lower iterations and smaller kappa.

q = waitbar(0, 'Filtering...');

% just for single set, later to be "for looped" for all sets    
j=1;  

%number of frames
numphases=size(Recon(j).M,3);

for i=1:numphases
         
        Recon(j).Vx_SM (:,:,i) = anisodiff2D(Recon(j).Vx(:,:,i),...
            Recon(j).M(:,:,i),num_iter,delta_t, kappa,option);
        Recon(j).Vy_SM (:,:,i) = anisodiff2D(Recon(j).Vy(:,:,i),...
            Recon(j).M(:,:,i),num_iter,delta_t, kappa,option);
        Recon(j).Vz_SM (:,:,i) = anisodiff2D(Recon(j).Vz(:,:,i),...
            Recon(j).M(:,:,i),num_iter,delta_t, kappa,option);

         waitbar(i/numphases,q,'Filtering...');      
end

save cs_processed.mat Recon
close(q)


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

