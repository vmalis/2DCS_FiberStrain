%=========================================================================
%   Step4 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strain Rate Toolkit
%=========================================================================
%
%
% INput: (output of Step2_ROIs)
%           1) results.mat
%
% OUTput:   1) Data_array.mat        multidimensional array of all data
%           2) Processed_data.mat    overal structure
%           3) Extracted_data.mat    at peak values and same level of force
%_____________________________________________________
% required subroutines:
%       1)rdir
%_____________________________________________________
% 
% v.2.0 
% (improved from v.1.0 used for Age-related Differences in Strain Rate Tensor paper (MRM 2014))
% 
% written by Vadim Malis
% 09/15 at UCSD RIL
%==========================================================================
%
% More info in HOWTO_2dSR.pdf

%% Read all the data

window_info=sprintf('Choose Data_VEPC folder');
PathName = uigetdir('~/Desktop',window_info);
cd(PathName)

file_list=rdir('**/results.mat');



for i=1:size(file_list,1)

    PathName = file_list(i).name;
    PathName = PathName(1:end-12);
    cd(PathName)   
    
    
    load('results.mat')
    load('force_data_calibrated.mat')

    DATA(i).patient_id=PatientID;
    DATA(i).series_name=Series_name;
    DATA(i).location=SliceLocation;
    DATA(i).sr_angle=ANGL;
    %DATA(i).fiber_angle=Fiber_ANGL;
    DATA(i).nev=NEV;
    DATA(i).pev=PEV;
    DATA(i).sev=SEV;
    DATA(i).quality=quality;
    DATA(i).force=force_mean;
    DATA(i).mvc=MAX_force;
    DATA(i).max_force=max(force_mean);
    DATA(i).frame_max_force=frame;
    DATA(i).force_level=0;
    DATA(i).force_level_frame=0;
    
    [~,temp]=max(DATA(1).force)
    DATA(i).frame_max_force_permuted=round(temp/600*22);
    
    for z=1:size(NEV,1)
    [~,index]=min(NEV(z,1:11));
    b(z)=index;
    end
    
    DATA(i).peak_ev_frame=b;
    
    DATA(i).force_at_peak=0;
    DATA(i).post_peak_same_force_frame=0;
    cd ..
    cd ..
    cd ..
    cd ..
    cd ..
    
end

clearvars -except DATA i

A=nestedSortStruct(DATA,{'initials','timing','location'});
num_sub=size(unique({A.initials}),2);
DATA=A;

% Overal 5D Array
%--------------------------------------------------------------------------
% building arrays combining information for all dimensions

% all the arrays below mantain the following size:
% 1 - number of subjects (variable)
% 2 - timing             (3, pre, post, post-rehab)
% 3 - number of slices   (8, if less then filled with NAN: sorted left to right)
% 4 - rois               (3 rois: sorted top to bottom)
% 5 - number of frames   (22, fixed for this study)

% replacing initials with indexes 
initials=unique({A.initials});
A=struct_replace(A,'initials',initials,1:size(initials,2));

NEV     = nan(num_sub,3,8,3,22);
PEV     = nan(num_sub,3,8,3,22);
SEV     = nan(num_sub,3,8,3,22);
SR_ANGL = nan(num_sub,3,8,3,22);
%FR_ANGL = nan(num_sub,3,8,3,22);

timing_old=A(1).timing;
slice=1;

for i=1:size(A,2)
          
    timing_new=A(i).timing;
    
    if timing_new==timing_old
    
      NEV(A(i).initials,A(i).timing,slice,:,:)      = A(i).nev;
      PEV(A(i).initials,A(i).timing,slice,:,:)      = A(i).pev;
      SEV(A(i).initials,A(i).timing,slice,:,:)      = A(i).sev;
      SR_ANGL(A(i).initials,A(i).timing,slice,:,:)  = A(i).sr_angle;
%       FR_ANGL(A(i).initials,A(i).timing,slice,:,:)  = A(i).fiber_angle;  
    
      slice=slice+1;
     
    else

      slice=1;
      
      NEV(A(i).initials,A(i).timing,slice,:,:)      = A(i).nev;
      PEV(A(i).initials,A(i).timing,slice,:,:)      = A(i).pev;
      SEV(A(i).initials,A(i).timing,slice,:,:)      = A(i).sev;
      SR_ANGL(A(i).initials,A(i).timing,slice,:,:)  = A(i).sr_angle;
%       FR_ANGL(A(i).initials,A(i).timing,slice,:,:)  = A(i).fiber_angle;
      
      slice=slice+1;
      timing_old=timing_new;  
      
    end
  
            
end
%--------------------------------------------------------------------------

save('Data_array.mat','NEV')
save('Data_array.mat','PEV','-append')
save('Data_array.mat','SEV','-append')
save('Data_array.mat','SR_ANGL','-append')
save('Data_array.mat','FR_ANGL','-append')

clearvars -except DATA


%--------------------------------------------------------------------------
%find max force of weakest attempt of each subject
force=[DATA.max_force];

[~,Unique_initials_index]=unique({DATA.initials});
Unique_initials_index=[Unique_initials_index;size(DATA,2)+1];


for i=1:size(Unique_initials_index,1)-1
    [DATA(Unique_initials_index(i):Unique_initials_index(i+1)-1).force_level]=deal(min(force(Unique_initials_index(i):Unique_initials_index(i+1)-1)));
end
    
%find frame corresponding to max weak force

for i=1:size(DATA,2)
  
[~,a]=min(abs(DATA(i).force(1:size(DATA(i).force,2)/2)-DATA(i).force_level));
frame=round(a/size(DATA(i).force,2)*size(DATA(i).nev,2));
DATA(i).force_level_frame=frame;
    
end



%for post find average frame of peak ev along the fiber

for i=1:size(DATA,2)
    
    if DATA(i).timing==2
    DATA(i).force_at_peak = [ DATA(i).force(round(DATA(i).peak_ev_frame(1)/22*600)), DATA(i).force(round(DATA(i).peak_ev_frame(2)/22*600)), DATA(i).force(round(DATA(i).peak_ev_frame(3)/22*600))];
    end

end


%finding number of slices in each case
% a unique patient IDs
% b number of slices in acquistion
% c index of a new series


a=unique({DATA.patient_id},'stable');
b=cellfun(@(x) sum(ismember({DATA.patient_id},x)),a,'un',0);
c=b;
c{1}=1;
for i=1:size(c,2)-1
c{i+1}=c{i}+b{i};
end

for i=1:size(c,2)
    if DATA(c{i}).force_at_peak==0
    INITIALS = DATA(c{i}).initials;
    
    bool=strfind({DATA.initials},INITIALS);
    index0 = find(~cellfun(@isempty,bool));

    index1_temp=find([DATA(index0).timing]==1);
    index1=index0(index1_temp);
    index2_temp=find([DATA(index0).timing]==2);
    index2=index0(index2_temp);
    index3_temp=find([DATA(index0).timing]==3);
    index3=index0(index3_temp);
    
    if size(index2,2)>size(index1,2)
       index2(size(index1,2)+1:end)=[];
    elseif size(index2,2)<size(index1,2)
       index1(size(index2,2)+1:end)=[];
    end    
     
    for j=1:size(index2,2)
    DATA(index1(j)).force_at_peak=DATA(index2(j)).force_at_peak;
    end
    
    if isempty(index3)==0
     
    if size(index3,2)>size(index1,2)
       index3(size(index1,2)+1:end)=[];
    elseif size(index3,2)<size(index1,2)
       index1(size(index3,2)+1:end)=[];
    end    
    
    for j=1:size(index3,2)
    DATA(index3(j)).force_at_peak=DATA(index2(j)).force_at_peak;
    end
    
    end
    
    end
end


for i=1:size(DATA,2)

    DATA(i).post_peak_same_force_frame=DATA(i).peak_ev_frame;
    
    
    if DATA(i).force_at_peak==0
    DATA(i).post_peak_same_force_frame=NaN;
    DATA(i).force_at_peak=NaN;
    end    
    
    if DATA(i).timing~=2 && size(DATA(i).force_at_peak,2)==3  
        for j=1:3    
            [~,temp]=min(abs(DATA(i).force(1:300)-DATA(i).force_at_peak(j)));
            temp=round(temp/600*22);
            if temp==0
                temp=1;
            end
            DATA(i).post_peak_same_force_frame(j)=temp;
        end
    end
      
        
        
end
save('Processed_data.mat','DATA');


%--------------------------------------------------------------------------
%% Create strucutre with data at peak eigen value along the fiber
DATA_peak(1).initials=0;
DATA_peak(1).timing=0;
DATA_peak(1).patient_id=0;
DATA_peak(1).location_sag=0;
DATA_peak(1).location_axial=0;
DATA_peak(1).sr_angle=0;
% DATA_peak(1).fiber_angle=0;
DATA_peak(1).nev=0;
DATA_peak(1).pev=0;
DATA_peak(1).sev=0;

for i=1:size(DATA,2)
    for j=1:size(DATA(i).nev,1)
    DATA_peak_temp(j).initials=DATA(i).initials;
    DATA_peak_temp(j).timing=DATA(i).timing; 
    DATA_peak_temp(j).patient_id=strcat(DATA(i).patient_id,num2str(j));
    DATA_peak_temp(j).location_sag=DATA(i).location;
    DATA_peak_temp(j).location_axial=j;
    DATA_peak_temp(j).sr_angle=DATA(i).sr_angle(j,DATA(i).peak_ev_frame(j));
%     DATA_peak_temp(j).fiber_angle=DATA(i).fiber_angle(j,DATA(i).peak_ev_frame(j));
    DATA_peak_temp(j).nev=DATA(i).nev(j,DATA(i).peak_ev_frame(j));
    DATA_peak_temp(j).pev=DATA(i).pev(j,DATA(i).peak_ev_frame(j));
    DATA_peak_temp(j).sev=DATA(i).sev(j,DATA(i).peak_ev_frame(j));
    end
    
    DATA_peak=[DATA_peak DATA_peak_temp];

end

DATA_peak(1)=[];

%--------------------------------------------------------------------------
%% Create strucutre with data at same force level
DATA_sfl(1).initials=0;
DATA_sfl(1).timing=0;
DATA_sfl(1).patient_id=0;
DATA_sfl(1).location_sag=0;
DATA_sfl(1).location_axial=0;
DATA_sfl(1).sr_angle=0;
DATA_sfl(1).fiber_angle=0;
DATA_sfl(1).nev=0;
DATA_sfl(1).pev=0;
DATA_sfl(1).sev=0;

for i=1:size(DATA,2)
    for j=1:size(DATA(i).nev,1)
    DATA_sfl_temp(j).initials=DATA(i).initials;
    DATA_sfl_temp(j).timing=DATA(i).timing; 
    DATA_sfl_temp(j).patient_id=strcat(DATA(i).patient_id,num2str(j));
    DATA_sfl_temp(j).location_sag=DATA(i).location;
    DATA_sfl_temp(j).location_axial=j;
    DATA_sfl_temp(j).sr_angle=DATA(i).sr_angle(j,DATA(i).force_level_frame);
    DATA_sfl_temp(j).fiber_angle=DATA(i).fiber_angle(j,DATA(i).force_level_frame);
    DATA_sfl_temp(j).nev=DATA(i).nev(j,DATA(i).force_level_frame);
    DATA_sfl_temp(j).pev=DATA(i).pev(j,DATA(i).force_level_frame);
    DATA_sfl_temp(j).sev=DATA(i).sev(j,DATA(i).force_level_frame);
    end
    
    DATA_sfl=[DATA_sfl DATA_sfl_temp];

end

DATA_sfl(1)=[];

clearvars -except DATA DATA_peak DATA_sfl


%--------------------------------------------------------------------------
%% Create strucutre with data at same force level as peak post
DATA_sfl_PP(1).initials=0;
DATA_sfl_PP(1).timing=0;
DATA_sfl_PP(1).patient_id=0;
DATA_sfl_PP(1).location_sag=0;
DATA_sfl_PP(1).location_axial=0;
DATA_sfl_PP(1).sr_angle=0;
DATA_sfl_PP(1).fiber_angle=0;
DATA_sfl_PP(1).nev=0;
DATA_sfl_PP(1).pev=0;
DATA_sfl_PP(1).sev=0;

for i=1:size(DATA,2)
    for j=1:size(DATA(i).nev,1)
    DATA_sfl_PP_temp(j).initials=DATA(i).initials;
    DATA_sfl_PP_temp(j).timing=DATA(i).timing; 
    DATA_sfl_PP_temp(j).patient_id=strcat(DATA(i).patient_id,num2str(j));
    DATA_sfl_PP_temp(j).location_sag=DATA(i).location;
    DATA_sfl_PP_temp(j).location_axial=j;
    DATA_sfl_PP_temp(j).sr_angle=DATA(i).sr_angle(j,DATA(i).post_peak_same_force_frame(j));
    DATA_sfl_PP_temp(j).fiber_angle=DATA(i).fiber_angle(j,DATA(i).post_peak_same_force_frame(j));
    DATA_sfl_PP_temp(j).nev=DATA(i).nev(j,DATA(i).post_peak_same_force_frame(j));
    DATA_sfl_PP_temp(j).pev=DATA(i).pev(j,DATA(i).post_peak_same_force_frame(j));
    DATA_sfl_PP_temp(j).sev=DATA(i).sev(j,DATA(i).post_peak_same_force_frame(j));
    end
    
    DATA_sfl_PP=[DATA_sfl_PP DATA_sfl_PP_temp];

end

DATA_sfl_PP(1)=[];

clearvars -except DATA DATA_peak DATA_sfl_PP DATA_sfl






%--------------------------------------------------------------------------
%% Average along sagittal slices using patient_id field
DATA_peak=nestedSortStruct(DATA_peak,'patient_id');
DATA_sfl=nestedSortStruct(DATA_sfl,'patient_id');
DATA_sfl_PP=nestedSortStruct(DATA_sfl_PP,'patient_id');


[~,Unique_patientID_index]=unique({DATA_peak.patient_id});
Unique_patientID_index=[Unique_patientID_index;size(DATA_peak,2)+1];


for i=1:size(Unique_patientID_index,1)-1

    DATA_peak_averaged(i).initials=DATA_peak(Unique_patientID_index(i)).initials;
    DATA_peak_averaged(i).timing=DATA_peak(Unique_patientID_index(i)).timing;
    DATA_peak_averaged(i).location_axial=DATA_peak(Unique_patientID_index(i)).location_axial;
    DATA_peak_averaged(i).sr_angle=mean([DATA_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).sr_angle]);
    DATA_peak_averaged(i).fiber_angle=mean([DATA_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).fiber_angle]);
    DATA_peak_averaged(i).nev=mean([DATA_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).nev]);
    DATA_peak_averaged(i).pev=mean([DATA_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).pev]);
    DATA_peak_averaged(i).sev=mean([DATA_peak(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).sev]);

%--------------------------------------------------------------------------

    DATA_sfl_averaged(i).initials=DATA_sfl(Unique_patientID_index(i)).initials;
    DATA_sfl_averaged(i).timing=DATA_sfl(Unique_patientID_index(i)).timing;
    DATA_sfl_averaged(i).location_axial=DATA_sfl(Unique_patientID_index(i)).location_axial;
    DATA_sfl_averaged(i).sr_angle=mean([DATA_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).sr_angle]);
    DATA_sfl_averaged(i).fiber_angle=mean([DATA_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).fiber_angle]);
    DATA_sfl_averaged(i).nev=mean([DATA_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).nev]);
    DATA_sfl_averaged(i).pev=mean([DATA_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).pev]);
    DATA_sfl_averaged(i).sev=mean([DATA_sfl(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).sev]);  
    
%--------------------------------------------------------------------------

    DATA_sfl_PP_averaged(i).initials=DATA_sfl_PP(Unique_patientID_index(i)).initials;
    DATA_sfl_PP_averaged(i).timing=DATA_sfl_PP(Unique_patientID_index(i)).timing;
    DATA_sfl_PP_averaged(i).location_axial=DATA_sfl_PP(Unique_patientID_index(i)).location_axial;
    DATA_sfl_PP_averaged(i).sr_angle=mean([DATA_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).sr_angle]);
    DATA_sfl_PP_averaged(i).fiber_angle=mean([DATA_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).fiber_angle]);
    DATA_sfl_PP_averaged(i).nev=mean([DATA_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).nev]);
    DATA_sfl_PP_averaged(i).pev=mean([DATA_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).pev]);
    DATA_sfl_PP_averaged(i).sev=mean([DATA_sfl_PP(Unique_patientID_index(i):Unique_patientID_index(i+1)-1).sev]);  
    
end

DATA_peak_averaged=nestedSortStruct(DATA_peak_averaged,'initials');
DATA_sfl_averaged=nestedSortStruct(DATA_sfl_averaged,'initials');
DATA_sfl_PP_averaged=nestedSortStruct(DATA_sfl_PP_averaged,'initials');

save('Extracted_data.mat','DATA_peak')
save('Extracted_data.mat','DATA_peak_averaged','-append')
save('Extracted_data.mat','DATA_sfl','-append')
save('Extracted_data.mat','DATA_sfl_averaged','-append')
save('Extracted_data.mat','DATA_sfl_PP','-append')
save('Extracted_data.mat','DATA_sfl_PP_averaged','-append')


