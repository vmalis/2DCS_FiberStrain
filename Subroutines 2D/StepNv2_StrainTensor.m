%=========================================================================
%   StepNv2 for 2d Strain Rate Analysis Toolbox (master script)
%
%       this script automaticaly does all aditional steps needed for
%       calucluation of strain rate tensor components in the basis of max
%       sheer component
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step2_ROIs)
%
%           1)  ev_calculated.mat
%           2)  mask.mat
%           4)  strain rate eigen values from ev_calculated_data.mat
%           6)  fiber data
%           7)* optional progress file allowing you to strat from the 
%               folder where you've ended up priveously
%
% OUTput:   1) modified result.mat (strain tensor in fiber basis added)
%
%_____________________________________________________
% required subroutines:
%       1)rdir
%_____________________________________________________
% 
% v.2.0 tailored for ULLS data analysis
% (improved from v.1.0 used for Age-related Differences in Strain Rate Tensor paper (MRM 2014))
% 
% written by Vadim Malis
% 05/16 at UCSD RIL
%==========================================================================
%
% More info in HOWTO_2dSR.pdf
%

%% load all set of possible data

window_info=sprintf('Choose Data_VEPC folder');
Path = uigetdir('~/Desktop',window_info);
cd(Path)
file_list_evs=rdir('**/ev_calculated_data.mat');
file_list_mask=rdir('**/mask.mat');
file_list_results=rdir('**/results.mat');


% loop (going through all data folders)

for i=1:size(file_list_results,1)
%for i=1
    
    % load mask
    load(file_list_mask(i).name)
    
    % load EV volumes
    load(file_list_evs(i).name,'EV_Negative')
    load(file_list_evs(i).name,'EV_Positive')
    load(file_list_evs(i).name, 'SR_Angle')
        
    % load fiber angle
    load(file_list_results(i).name,'Fiber_ANGL')
    load(file_list_results(i).name,'PatientID')
    load(file_list_results(i).name,'Series_name')
    load(file_list_results(i).name,'SliceLocation')
    
    % mask EVs
    
    SR_tensor=zeros(2,2,22,3);
    
    
    teta=deg2rad(45);
    
    for j=1:size(MASK,4)
      
         mask=squeeze(MASK(:,:,:,j));
         NEV=EV_Negative.*mask;
         PEV=EV_Positive.*mask;
         ANGLE=SR_Angle.*mask;
         
         [i_1,i_2,i_3] = ind2sub(size(mask),find(mask));
         
         SR_temp=zeros(2,2,size(i_3,1));
         
         for k=1:size(i_1,1)       
             %teta=Fiber_ANGL(j,i_3(k))-ANGLE(i_1(k),i_2(k),i_3(k));
             %teta=deg2rad(teta);
             TransformationMatrix=[cos(teta), -sin(teta);sin(teta), cos(teta)];
             SR=[PEV(i_1(k),i_2(k),i_3(k)),0;0,NEV(i_1(k),i_2(k),i_3(k))];
             SR_transformed=TransformationMatrix*SR*TransformationMatrix';
             SR_temp(:,:,k)=SR_transformed;
         end   
             
         [a,b]=unique(i_3);
         b=[b;size(i_3,1)];
         
         for k=1:size(a,1)
             SR_tensor(:,:,k,j)=mean(SR_temp(:,:,b(k):b(k+1)),3);
         end
         
    end
     
    
    PathName = file_list_results(i).name;
    PathName = PathName(1:end-12);
    cd(PathName)   
    
    if strcmp(PathName(9),'-')
        SR_data(i).initials= PathName(16:17);
        SR_data(i).timing  = 3;
    elseif strcmp(PathName(8),'/')
        SR_data(i).initials= PathName(9:10);
        SR_data(i).timing  = 1;
    elseif strcmp(PathName(9),'/')
        SR_data(i).initials= PathName(10:11);
        SR_data(i).timing  = 2;
    end
    
    
    
    
    SR_data(i).patient_id=PatientID;
    SR_data(i).series_name=Series_name;
    SR_data(i).location=SliceLocation;
    SR_data(i).SR_tensor=SR_tensor(:,:,:,:);
   
cd(Path)    
    
% end loop    
end
    
% create strain tensor results file (similar to results)
save('SR_tensor_max_sheer.mat','SR_data');    

load('Processed_data.mat')
SR_data_sorted=nestedSortStruct(SR_data,{'initials','timing','location'});

%frame information is used from the DATA structure
%extracting frame specific information and creating two new structures


%%

%--------------------------------------------------------------------------
%% Create strucutre with data at peak eigen value along the fiber
SR_peak(1).initials=0;
SR_peak(1).timing=0;
SR_peak(1).patient_id=0;
SR_peak(1).location_sag=0;
SR_peak(1).location_axial=0;
SR_peak(1).SR_ff=0;
SR_peak(1).SR_cc=0;
SR_peak(1).SR_fc=0;


for i=1:size(DATA,2)
    for j=1:size(DATA(i).nev,1)
    SR_peak_temp(j).initials=DATA(i).initials;
    SR_peak_temp(j).timing=DATA(i).timing; 
    SR_peak_temp(j).patient_id=strcat(DATA(i).patient_id,num2str(j));
    SR_peak_temp(j).location_sag=DATA(i).location;
    SR_peak_temp(j).location_axial=j;
    SR_peak_temp(j).SR_ff=SR_data_sorted(i).SR_tensor(1,1,DATA(i).peak_ev_frame(j),j);
    SR_peak_temp(j).SR_cc=SR_data_sorted(i).SR_tensor(2,2,DATA(i).peak_ev_frame(j),j);
    SR_peak_temp(j).SR_fc=SR_data_sorted(i).SR_tensor(1,2,DATA(i).peak_ev_frame(j),j);
    end
    
    SR_peak=[SR_peak SR_peak_temp];

end

SR_peak(1)=[];

%--------------------------------------------------------------------------
%% Create strucutre with data at same force level
SR_sfl(1).initials=0;
SR_sfl(1).timing=0;
SR_sfl(1).patient_id=0;
SR_sfl(1).location_sag=0;
SR_sfl(1).location_axial=0;
SR_sfl(1).SR_ff=0;
SR_sfl(1).SR_cc=0;
SR_sfl(1).SR_fc=0;


for i=1:size(DATA,2)
    for j=1:size(DATA(i).nev,1)
    SR_sfl_temp(j).initials=DATA(i).initials;
    SR_sfl_temp(j).timing=DATA(i).timing; 
    SR_sfl_temp(j).patient_id=strcat(DATA(i).patient_id,num2str(j));
    SR_sfl_temp(j).location_sag=DATA(i).location;
    SR_sfl_temp(j).location_axial=j;
    SR_sfl_temp(j).SR_ff=SR_data_sorted(i).SR_tensor(1,1,DATA(i).force_level_frame,j);
    SR_sfl_temp(j).SR_cc=SR_data_sorted(i).SR_tensor(2,2,DATA(i).force_level_frame,j);
    SR_sfl_temp(j).SR_fc=SR_data_sorted(i).SR_tensor(1,2,DATA(i).force_level_frame,j);
    end
    
    SR_sfl=[SR_sfl SR_sfl_temp];

end

SR_sfl(1)=[];

clearvars -except DATA SR_peak SR_sfl SR_data_sorted


%--------------------------------------------------------------------------
%% Create strucutre with data at same force level as peak post
SR_sfl_PP(1).initials=0;
SR_sfl_PP(1).timing=0;
SR_sfl_PP(1).patient_id=0;
SR_sfl_PP(1).location_sag=0;
SR_sfl_PP(1).location_axial=0;
SR_sfl_PP(1).SR_ff=0;
SR_sfl_PP(1).SR_cc=0;
SR_sfl_PP(1).SR_fc=0;

for i=1:size(DATA,2)
    for j=1:size(DATA(i).nev,1)
    SR_sfl_PP_temp(j).initials=DATA(i).initials;
    SR_sfl_PP_temp(j).timing=DATA(i).timing; 
    SR_sfl_PP_temp(j).patient_id=strcat(DATA(i).patient_id,num2str(j));
    SR_sfl_PP_temp(j).location_sag=DATA(i).location;
    SR_sfl_PP_temp(j).location_axial=j;
    SR_sfl_PP_temp(j).SR_ff=SR_data_sorted(i).SR_tensor(1,1,DATA(i).post_peak_same_force_frame(j),j);
    SR_sfl_PP_temp(j).SR_cc=SR_data_sorted(i).SR_tensor(2,2,DATA(i).post_peak_same_force_frame(j),j);
    SR_sfl_PP_temp(j).SR_fc=SR_data_sorted(i).SR_tensor(1,2,DATA(i).post_peak_same_force_frame(j),j);
    end
    
    SR_sfl_PP=[SR_sfl_PP SR_sfl_PP_temp];

end

SR_sfl_PP(1)=[];

clearvars -except DATA SR_peak SR_sfl_PP SR_sfl SR_data_sorted


%--------------------------------------------------------------------------
%% Average along sagittal slices using patient_id field
SR_peak=nestedSortStruct(SR_peak,'patient_id');
SR_sfl=nestedSortStruct(SR_sfl,'patient_id');
SR_sfl_PP=nestedSortStruct(SR_sfl_PP,'patient_id');


[~,Unique_patientID_index]=unique({SR_peak.patient_id});
Unique_patientID_index=[Unique_patientID_index;size(SR_peak,2)+1];


for i=1:size(Unique_patientID_index,1)-1

    SR_peak_averaged(i).initials=SR_peak(Unique_patientID_index(i)).initials;
    SR_peak_averaged(i).timing=SR_peak(Unique_patientID_index(i)).timing;
    SR_peak_averaged(i).location_axial=SR_peak(Unique_patientID_index(i)).location_axial;
    SR_peak_averaged(i).SR_ff=stdfilter([SR_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_ff]);
    SR_peak_averaged(i).SR_cc=stdfilter([SR_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_cc]);
    SR_peak_averaged(i).SR_fc=stdfilter([SR_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_fc]);    
    
%--------------------------------------------------------------------------

    SR_sfl_averaged(i).initials=SR_sfl(Unique_patientID_index(i)).initials;
    SR_sfl_averaged(i).timing=SR_sfl(Unique_patientID_index(i)).timing;
    SR_sfl_averaged(i).location_axial=SR_sfl(Unique_patientID_index(i)).location_axial;
    SR_sfl_averaged(i).SR_ff=stdfilter([SR_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_ff]);
    SR_sfl_averaged(i).SR_cc=stdfilter([SR_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_cc]);
    SR_sfl_averaged(i).SR_fc=stdfilter([SR_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_fc]);
    
%--------------------------------------------------------------------------

    SR_sfl_PP_averaged(i).initials=SR_sfl_PP(Unique_patientID_index(i)).initials;
    SR_sfl_PP_averaged(i).timing=SR_sfl_PP(Unique_patientID_index(i)).timing;
    SR_sfl_PP_averaged(i).location_axial=SR_sfl_PP(Unique_patientID_index(i)).location_axial;
    SR_sfl_PP_averaged(i).SR_ff=stdfilter([SR_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_ff]);
    SR_sfl_PP_averaged(i).SR_cc=stdfilter([SR_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_cc]);
    SR_sfl_PP_averaged(i).SR_fc=stdfilter([SR_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).SR_fc]);
    
end

SR_peak_averaged=nestedSortStruct(SR_peak_averaged,'initials');
SR_sfl_averaged=nestedSortStruct(SR_sfl_averaged,'initials');
SR_sfl_PP_averaged=nestedSortStruct(SR_sfl_PP_averaged,'initials');

save('Extracted_max_sheer_SR_data.mat','SR_peak')
save('Extracted_max_sheer_SR_data.mat','SR_peak_averaged','-append')
save('Extracted_max_sheer_SR_data.mat','SR_sfl','-append')
save('Extracted_max_sheer_SR_data.mat','SR_sfl_averaged','-append')
save('Extracted_max_sheer_SR_data.mat','SR_sfl_PP','-append')
save('Extracted_max_sheer_SR_data.mat','SR_sfl_PP_averaged','-append')
%finaly preparing arrays for plots



