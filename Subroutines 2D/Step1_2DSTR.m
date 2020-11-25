%=========================================================================
%   Step1 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strane Rate Toolkit
%=========================================================================
%
% This script
%     a) sorts mr images
%     b) applies eddy current correction (subtraction off average image)
%     c) performs (2dAnisoDif) filtering
%     d) reads and saves header
%     e) calculates and saves avarage force
%     f) performs permutation for max force to be in the middle of cycle
%     g) eigenvalue and SR angle maps
%     h) reads and sorts fiber endpoints information
%
% INput: *.dcm    dicom series from the GE MR scanner
%        *.csv    force recording
%        *.csv    mvc force recording
%        *.dcm    highres FGRE used to identify fiber endpoints
%        *.csv    fiber endpoints
%
% OUTput: 1) sorted velocities stack        |
%         2) filtered velocities stack      |->  mri_sorted.mat
%         3) header information             |
%         4) analysed force data            -    force_data.mat; graphs
%         5) eigen values stacks            ?    positive_eig.dat; negative_eig.dat; sum_eig.dat   
%         6) SR angle stack                 ?    angle.dat
%         7) data for plot                  ?    ev_calculated_data.mat
%__________________________________________________________________________
% required subroutines:
%
%   1)vepc_2d_imsort
%   2)force_data
%   3)im_subtract
%   4)im_mean
%   5)anisodiff2D
%   6)multiWaitbar
%   7)fiber_read
%   8)fiber_end
%
%__________________________________________________________________________
% 
% v.2.0 
% (improved from v.1.0 used for Age-Related Differences in Strain Rate Tensor paper (MRM 2014) )
% 
% written by Vadim Malis
% 09/14 at UCSD RIL
%==========================================================================


% Select the folder where all the folders with the subjects are
% folder structure should be according to Step0_FStructure.m

window_info=sprintf('pick root data folder');
PathName = uigetdir('~/Desktop',window_info);
cd(PathName)

% folder_list function is valid for MAC only

[path_study,study_n]=folder_list(pwd);


multiWaitbar('Overal progress...', 0, 'Color', 'r');

for counter_study=1:study_n

cd(path_study(counter_study).name)

[path_subject,subject_n]=folder_list(pwd);


multiWaitbar('per subject...', 0, 'Color', 'y');

for counter_subject=1:subject_n

cd(path_subject(counter_subject).name)
% cd('patient')
% cd('study')
[path_series,k]=folder_list(pwd);    


multiWaitbar('per series...', 0, 'Color', 'b');

for counter1=1:k                %per series

cd(path_series(counter1).name)
dicom_path=dir('*dcm');    

cd


%% Read original images and sort


if counter1==1
    
   I=mat2gray(dicomread(dicom_path(counter1).name));
   figure
   imshow(I)
    
    choice = questdlg('Is th orientation left or right sided?', ...
	'Orientation', ...
	'Left','Right','Right');
     % Handle response
        switch choice
            case 'Left'            
                orientation=1;        
            case 'Right'
                orientation=0;
                
        end
    close
end

im_data=vepc_2d_imsort(dicom_path,orientation);
load(im_data);

%% ------------------------------------------------------------------------
% Force analysis
bpm=20;
path=pwd;
force_filename=force_analysis(path,numphases,bpm,1);
load(force_filename);


%% ------------------------------------------------------------------------
% Calculating SR and angle maps

% Magnitude image file 
magimage1 = im_m(:,:,1);

%   Allocate the gradient image buffers
%   Create 4 matrices to store the EigenValues of Every Voxel
f1 = zeros(size(v_ap),'single');
f2 = zeros(size(v_ap),'single');
f3 = zeros(size(v_ap),'single');
f4 = zeros(size(v_ap),'single');

%   Allocate arrays for filtered velocities
v_rl_sm=zeros(size(v_rl));
v_ap_sm=zeros(size(v_rl));
v_si_sm=zeros(size(v_rl));

%---------2D AnisotropicDiffusion Filter-----------------------------------

        num_iter=10; kappa=2; option=1; delta_t=1/7;
        % calling the anisodiff2d code; options are set above; 
        % less smoothing with lower iterations and smaller kappa.
       
        multiWaitbar('Filtering...', 0, 'Color', 'g');
        
for i=1:numphases
         
        v_rl_sm(:,:,i) = anisodiff2D(v_rl(:,:,i),...
            im_m(:,:,i),num_iter,delta_t, kappa,option);
        v_ap_sm(:,:,i) = anisodiff2D(v_ap(:,:,i),...
            im_m(:,:,i),num_iter,delta_t, kappa,option);
        v_si_sm(:,:,i) = anisodiff2D(v_si(:,:,i),...
            im_m(:,:,i),num_iter,delta_t, kappa,option);
        
%----------Spatial derivative----------------------------------------------
%   2D_Strain Tensor    
         [f1(:,:,i), f2(:,:,i)] = gradient(v_ap_sm(:,:,i));
         [f3(:,:,i), f4(:,:,i)] = gradient(v_si_sm(:,:,i));
         multiWaitbar('Filtering...', i/numphases);      
end
     
multiWaitbar('Filtering...', 'Close');  

%% ------------------------------------------------------------------------
% Prealocation of variables for Eigen Value decomposition

% Create a Matrix to store the 4 components of the StrainTensor of every
% [STxx,STxy,STyx, STyy]
StrainT = zeros([size(v_ap) 4],'single');

% Create 2 matrices to store the EigenValues of every Voxel
Y_1 = zeros(size(v_ap_sm),'single');
Y_2 = zeros(size(v_ap_sm) ,'single');
Y_sqrt = zeros(size(v_ap_sm) ,'single');
Y_sum = zeros(size(v_ap_sm) ,'single');

% Create a maxtrix to store the (main) strain direction in each pixel
VectorF_red=zeros([size(v_ap_sm) 2],'single');
VectorF_blue=zeros([size(v_ap_sm) 2],'single');

%
angle=zeros(size(v_ap_sm),'single');
x1angle=zeros(size(v_ap_sm),'single');
y1angle=zeros(size(v_ap_sm),'single');
x2angle=zeros(size(v_ap_sm),'single');
y2angle=zeros(size(v_ap_sm),'single');
angle_1=zeros(size(v_ap_sm),'single');

%% ------------------------------------------------------------------------
% Eigenvalue problem CALCS

multiWaitbar('EV decomposition...', 0, 'Color', 'g');


A=size(v_ap_sm,1);
B=size(v_ap_sm,2);
C=size(v_ap_sm,3);

for a=1:A
    for b=1:B
        for c=1:C
            
%             Calculating the Strain Tensor, Transpose, Eigenvectors and
%             values at each voxel, Sorting and Ordering Eigenvalues from
%             samllest to largest (Normalizing each gradient by its voxel
%             dimensions as well)

            if (magimage1(a,b))>0.01
                
            StrainTensor=[(f1(a,b,c))*(resolution) (f2(a,b,c))*(resolution);(f3(a,b,c))*(resolution) (f4(a,b,c))*(resolution)];
            StrainTensorTrans=[(f1(a,b,c))*(resolution) (f3(a,b,c))*(resolution);(f2(a,b,c))*(resolution) (f4(a,b,c))*(resolution)];
            L=0.5*(StrainTensor+StrainTensorTrans);
            [EigenVectors,S]=eig(L);EigenValues=diag(S);
            [t,index]=sort(EigenValues);
            
%Always second EigenValue is positive
            EigenValues=EigenValues(index);            
            EigenVectors = EigenVectors(:,index);   
            Y_1(a,b,c)=EigenValues(1)*1000;
            Y_2(a,b,c)=EigenValues(2)*1000;
            Y_sqrt(a,b,c)=sqrt(power(EigenValues(1),2)+power(EigenValues(2),2))*1000;
            Y_sum(a,b,c) = Y_1(a,b,c)+ Y_2(a,b,c);
       
%Bring all EigenVectors to first and second quadrant

            %NEV
            if((EigenVectors(1,1)>0 && EigenVectors(2,1)>0)||(EigenVectors(1,1)<0 && EigenVectors(2,1)>0))
            VectorF_red(a,b,c,:)= EigenVectors(:,1);
            x1angle(a,b,c) = (180/pi)*acos(EigenVectors(1,1));
            elseif((EigenVectors(1,1)>0&&EigenVectors(2,1)<0)||(EigenVectors(1,1)<0 && EigenVectors(2,1)<0))
            VectorF_red(a,b,c,:)= -EigenVectors(:,1);
            x1angle(a,b,c) = 180-((180/pi)*acos(EigenVectors(1,1)));
            else
            end


            %PEV
            if((EigenVectors(1,2)>0 && EigenVectors(2,2)>0)||(EigenVectors(1,2)<0 && EigenVectors(2,2)>0))
            VectorF_blue(a,b,c,:)= EigenVectors(:,2);
            elseif((EigenVectors(1,2)>0 && EigenVectors(2,2)<0)||(EigenVectors(1,2)<0 && EigenVectors(2,2)<0))
            VectorF_blue(a,b,c,:)= -EigenVectors(:,2);
            else
            end
   
            end
        end
    end
    multiWaitbar('EV decomposition...', a/A);
    
end

multiWaitbar('EV decomposition...', 'Close');

%Angle to 1st and 4th quadrant (NEV with +x)
x1angle(x1angle>90)=(x1angle(x1angle>90)-180);
x1angle=(-1)*x1angle;


EV_Negative=Y_1;
EV_Positive=Y_2;
EV_Sum=Y_sum;



%% ------------------------------------------------------------------------
% Saving calculated volumes: .dat for convinience to preview in ImageJ
% For processing necessary data are saved in ev_calculated.mat

fid = fopen('negative_eig.dat','w+','l');
Y_1 = permute(Y_1,[2 1 3]);
fwrite(fid,Y_1,'float');
fclose(fid);
% 
fid = fopen('positive_eig.dat','w+','l');
Y_2 = permute(Y_2,[2 1 3]);
fwrite(fid,Y_2,'float');
fclose(fid);
% 
fid = fopen('sum_eig.dat','w+','l');
Y_sum = permute(Y_sum,[2 1 3]);
fwrite(fid,Y_sum,'float');
fclose(fid);
% 
fid = fopen('angle.dat','w+','l');
ANGLE = permute(x1angle,[2 1 3]);
fwrite(fid,ANGLE,'float');
fclose(fid);

% %save all workspace data
% save ('ev_calculated_data.mat')


SR_Angle=x1angle;
NEVector=VectorF_red;
PEVector=VectorF_blue;

%--------------------------------EV decomposition data---------------------
save ('ev_calculated_data.mat','EV_Negative');
save ('ev_calculated_data.mat','EV_Positive','-append');
save ('ev_calculated_data.mat','EV_Sum','-append');
save ('ev_calculated_data.mat','SR_Angle','-append');
save ('ev_calculated_data.mat','NEVector','-append')
save ('ev_calculated_data.mat','PEVector','-append')

%-------------------------------velocities---------------------------------
% save ('ev_calculated_data.mat','v_ap','-append');
% save ('ev_calculated_data.mat','v_rl','-append')
% save ('ev_calculated_data.mat','v_si','-append')
save ('ev_calculated_data.mat','v_ap_sm','-append')
save ('ev_calculated_data.mat','v_rl_sm','-append')
save ('ev_calculated_data.mat','v_si_sm','-append')

%---------------------------------header-----------------------------------
save ('ev_calculated_data.mat','PatientID','-append');
save ('ev_calculated_data.mat','Series_name','-append');
save ('ev_calculated_data.mat','RotationMatrix','-append');
save ('ev_calculated_data.mat','SliceLocation','-append');
save ('ev_calculated_data.mat','c','-append');
save ('ev_calculated_data.mat','r','-append');
save ('ev_calculated_data.mat','resolution','-append');
save ('ev_calculated_data.mat','dt','-append');
save ('ev_calculated_data.mat','numphases','-append');



cd ..

multiWaitbar('per series...', counter1/k);

end

cd ..
% cd .. %out from osirix's folder study
% cd .. %out from osirix's folder patient

multiWaitbar('per series...', 'Close');
multiWaitbar('per subject...', counter_subject/subject_n);

end

cd .. %out from subject initials folder

multiWaitbar('per subject...', 'Close');
multiWaitbar('Overal progress...', counter_study/study_n);

end

multiWaitbar('Overal progress...', 'Close');
clear all
clc


cd .. %out of study folder
%%cd .. %out of modality folder (ISO, ECC, ...)

pwd
%% ------------------------------------------------------------------------
% Fiber processing (obtain coordinates, save data)
if exist('DATA_FIBER')
    cd ('DATA_FIBER')
else

if exist(fullfile(cd, 'fiber data.mat'), 'file') == 2
 
else
    h = msgbox('Last step: read fiber information');
    fiber_read(pwd);
    clear all
    close
end
end