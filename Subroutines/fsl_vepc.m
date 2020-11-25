function [DWI,baseline_corrected,Grad,Location,info]=fsl_vepc
% part of 3dSR
%
% subroutine to run FSL eddy current and fieldmap correction
% runs eddy_fieldmap.sh
% _____________________________________________________
% written by Vadim Malis
% 10/17 at UCSD RIL


b=userpath;
%a=findstr(t,':');
%b=t(1:a(1)-1);

path1 = getenv('PATH');
path1 = [path1 ':' b '/Users/vadmalis/Documents/MATLAB/DTI_v2:/usr/local/fsl/bin'];
setenv('PATH', path1);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
datapath=uigetdir();
cd(datapath)

all_files = dir;
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir);

cd EPI
% read 1st image dicom series
dicom_filenames=dir('*.dcm');
info = dicominfo(dicom_filenames(1).name);
ID = info.PatientID;

% Read Gradients of 1st image
GradientDirectionA=[info.Private_0019_10bb, info.Private_0019_10bc, ...
     info.Private_0019_10bd];

% number of images in the DWI set
n=length(dicom_filenames);

% allocate memory for Gradient table
GradTable=zeros(n,3);
Location=zeros(n,1);

% waitbar
h.waitbar = waitbar(0,'Please be patient!');


% initial number of gradient direction configuration
i=0;


% processing all the files of the set
for k=2:n
   
    
    %fname=sprintf('%s%.4d.dcm',general_filename,k);
    fname=dicom_filenames(k).name;

    % read header
    info = dicominfo(fname);

    % obtaining gradient direction and location
    
    GradientDirectionB=[info.Private_0019_10bb, info.Private_0019_10bc,...
	info.Private_0019_10bd];
    Location(k)=info.SliceLocation;

    % writing gradient direction to the GradTable array
    GradTable(k,:) = GradientDirectionB;
    
    %waitbar
    waitbar(k/n)
    
    
    % Check if gradients change
    if GradientDirectionA == GradientDirectionB
       dummy=0;
       
    else
        i=i+1;
        GradientDirectionA = GradientDirectionB;
        
    end
    
end
close all

% grad table aray disposes duplicate lines
GradTable=unique(GradTable,'rows','stable');
Location(1)=[];
Location=unique(Location);

%Check phase encode direction and perform permutation if needed
if info.InPlanePhaseEncodingDirection=='ROW'
index_PhaseEncode_ROW=[2,1,3];
Grad=GradTable(:,index_PhaseEncode_ROW);
end

clearvars -except Grad Location ID num_dir
cd ..

%%
b=userpath;
%a=findstr(t,':');
%b=t(1:a(1)-1);

path1 = getenv('PATH');
path1 = [path1 ':' b '/DTI_v2:/FSL/fsl/bin'];
setenv('PATH', path1);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');


    %% call for c-shell script and perform eddy current
    !eddy_fieldmap_NFM.sh
    fname=dir('*.bval');
    fileID = fopen(fname.name,'r');
    formatSpec = '%f';
    b = fscanf(fileID,formatSpec);
    fclose(fileID);
    b(3:end)=[];
    b(1)=[];

    %% read header
    Info=nii_read_header('EPI.nii');
    nslices=Info.Dimensions(3);
    info.b=b;
    info.x=Info.QoffsetX;
    info.y=Info.QoffsetY;
    info.ID=ID;

    % read DWI corrected volume
    gunzip('DWI_EC.nii.gz')
    delete DWI_EC.nii.gz
    dwi_header=nii_read_header('DWI_EC.nii');
    DWI=nii_read_volume(dwi_header);
    DWI=flip(permute(DWI,[2,1,3,4]),3);
    DWI=double(flip(DWI,1));
    baseline_corrected=squeeze(DWI(:,:,:,1));
    baseline_corrected_temp = permute(reshape(baseline_corrected,[size(baseline_corrected),1]),[1,2,4,3]);

    % read original DWI
    epi_header=nii_read_header('EPI.nii');
    EPI=nii_read_volume(epi_header);
    EPI=flip(permute(EPI,[2,1,3,4]),3);
    baseline=squeeze(EPI(:,:,:,1));
    baseline=flip(baseline,1);
    baseline_temp =double(permute(reshape(baseline,[size(baseline),1]),[1,2,4,3]));
    clear EPI

    % %montage
    for i=1:nslices
        baseline_corrected_temp(:,:,:,i) = mat2gray(baseline_corrected_temp(:,:,:,i));
        baseline_temp(:,:,:,i) = mat2gray(baseline_temp(:,:,:,i));
    end

    %----1 original baseline
    montage(baseline_temp)
    export_fig('Montage_baseline.png', '-png','-m16')
    close

    %----2 baseline corrected
    montage(baseline_corrected_temp)
    export_fig('Montage_baseline_corrected.png', '-png','-m16')
    close
end
    
