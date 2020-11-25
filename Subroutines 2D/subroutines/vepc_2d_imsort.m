function [filename] = vepc_2d_imsort(path,orientation)

%==========================================================================
%  Subroutine to perform dicom image sort
%==========================================================================
%
% Read original images and sort (based on David Shin rdgeimst_all.m)
% 
% INput:        structure with list of *.dcm images in a folder
% OUTput:       *.mat files with sorted images and image parameters
%--------------------------------------------------------------------------
%
%   can be called directly:
%
%           PathName = uigetdir('~/Desktop');
%           cd(PathName)
%           dicom_path=dir('*dcm');
%           vepc_2d_imsort(dicom_path)
%
%--------------------------------------------------------------------------
% written by Vadim Malis (based on David Shin rdgimst_all.m)
% 09/14 at UCSD RIL
%==========================================================================



hdinfo = dicominfo(path(1).name);
resolution = 1 / (hdinfo.PixelSpacing(1) / 10);  %pixels/cm
numim = hdinfo.ImagesInAcquisition;
numphases = numim/5;
r = hdinfo.Rows;
c = hdinfo.Columns;
Series_name=hdinfo.SeriesDescription;
PatientID=hdinfo.PatientID;
SliceLocation=hdinfo.SliceLocation;
RotationMatrix_temp=reshape(hdinfo.ImageOrientationPatient,3,2);
RotationMatrix=cat(2,RotationMatrix_temp,cross(RotationMatrix_temp(:,1),...
    RotationMatrix_temp(:,2)));



triggertime = zeros(numim,1);
triggertime_std = zeros(numim,1);
venc = zeros(numim,1);
venc_axis = zeros(numim,1);

insnum = zeros(numim,1); %instance number
im = zeros(r,c,numim);
im_al = zeros(r,c,numphases);
im_m = zeros(r,c,numphases); 
im_rl = zeros(r,c,numphases);
im_ap = zeros(r,c,numphases);
im_si = zeros(r,c,numphases);


%% Loading images
multiWaitbar('Reading images', 0,'Color', 'g');

for i = 1:numim
    im(:,:,i) = double(dicomread(path(i).name));
    dcmhd = dicominfo(path(i).name);
    triggertime(i) = dcmhd.TriggerTime;
    venc_axis(i) = dcmhd.(dicomlookup('0019','10cb'));
    venc(i) = dcmhd.(dicomlookup('0019','10cc'));
    insnum(i) = dcmhd.InstanceNumber;
multiWaitbar('Reading images', i/numim);
end
multiWaitbar('Reading images', 'Close');


%% Sorting images

multiWaitbar('Sorting images', 0,'Color', 'g');

[insnum_std,index] = sort(insnum);

c_al = 1;
c_m = 1;
c_rl = 1;
c_ap = 1;
c_si = 1;

for i=1:numim
    
    if (i <= numphases)
         %velocity converted to cm/s & SI -> Positive
        im_al(:,:,c_al) = im(:,:,index(i)).*0.1; 
        c_al = c_al + 1;
    elseif (i > numphases) && (i <= numphases*2) 
        im_m(:,:,c_m) = im(:,:,index(i));
        c_m = c_m + 1;
    elseif (i > numphases*2) && (i <= numphases*3)
        im_rl(:,:,c_rl) = im(:,:,index(i)).*0.1; 
        c_rl = c_rl + 1;
    elseif (i > numphases*3) && (i <= numphases*4)
        im_ap(:,:,c_ap) = im(:,:,index(i)).*0.1; 
        c_ap = c_ap + 1;        
    else
        im_si(:,:,c_si) = im(:,:,index(i)).*0.1; 
        c_si = c_si + 1;        
    end 
    
triggertime_std(i) = triggertime(index(i));

% eddy current correction

error_rl = im_mean(im_rl);
error_ap = im_mean(im_ap);
error_si = im_mean(im_si);

v_rl = im_subtract(im_rl, error_rl);
v_ap = im_subtract(im_ap, error_ap);
v_si = im_subtract(im_si, error_si);


multiWaitbar('Sorting images', i/numim);


end


% find dt
dt = zeros(1,numphases-1);
for ii=2:numphases
    dt(ii-1) = (triggertime_std(ii) - triggertime_std(ii-1))/1000;
end



if orientation==1
    
    im_al=flip(im_al,2);
	im_ap=flip(im_ap,2);
	im_m=flip(im_m,2);
	im_rl=flip(im_rl,2);
	im_si=flip(im_si,2);
	v_ap=flip(v_ap,2);
	v_rl=flip(v_rl,2);
	v_si=flip(v_si,2);
                
	v_rl=v_rl*(-1);
	im_rl=im_rl*(-1);
                
end


multiWaitbar('Sorting images', 'Close');  
filename='mri_data.mat';

%% save variables as mat file
save (filename,'im_al')
save (filename,'im_m','-append')
save (filename,'im_rl','-append')
save (filename,'im_ap','-append')
save (filename,'im_si','-append')
save (filename,'v_rl','-append')
save (filename,'v_ap','-append')
save (filename,'v_si','-append')
save (filename,'insnum_std','-append')
save (filename,'dt','-append')
save (filename,'numphases','-append')
save (filename,'resolution','-append')
save (filename,'r','-append')
save (filename,'c','-append')
save (filename,'Series_name','-append')
save (filename,'PatientID','-append')
save (filename,'SliceLocation','-append')
save (filename,'RotationMatrix','-append')

end