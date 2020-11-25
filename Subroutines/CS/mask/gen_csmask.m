function [K_pattern,mask,pdf]=gen_csmask(nphases, nreadout, nframes,cs_factor)

% Generates undersampling pattern
%   
% Input:   
%      phase encoding size
%      readout size
%      number of frames
%      undersampling factor
%
% Output:  
%      K_pattern    - 1d pattern
%      CS_pattern   - mask
%      pdf          ? probability density function
%
% example call:
%  [K_pattern,CS_pattern,PDF]=gen_csmask(106,256,22,2);
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis

%% create undersampling pattern
phase_size=nphases;
readout_size=nreadout;
frames=nframes;
cs_factor=1/cs_factor;

file_name = strcat('CS_pattern','_x',num2str(1/cs_factor),'.txt');
file_name_mat = strcat('CS','_x',num2str(1/cs_factor),'.mat');
file_name_img = strcat('CS','_x',num2str(1/cs_factor),'.eps');

% deafualt polynomial power suggested from Lustig et al. n=20
poly_power=20;

% gen probability density function
[pdf,~]=genPDF([phase_size,1],poly_power,cs_factor,2,0,0);
K_pattern=zeros(frames,phase_size);
N=zeros(frames,1);
b=zeros(nframes,ceil(phase_size*cs_factor));

for i=1:frames

    [mask,~,n] = genSampling(pdf,3,1); 
    N(i)=n;
    a=find(mask)-1;  
        if size(a,2)>ceil(phase_size*cs_factor)
            mask(a(end)) = 0;
            a(end)=[];
        elseif size(a,2)<ceil(phase_size*cs_factor)
            
            for r=0:nphases
                qq=find(a==r, 1);
                if isempty(qq)
                    a=[a, r];
                break    
                end
            end
            
        end
    a=sort(a);    
  
    b(i,:)=a;
    K_pattern(i,:) = mask;

end

b=b';
b=b(:);
b=b';
dlmwrite(file_name,b,'-append','delimiter','\n')

CS_pattern=repmat(K_pattern,[1,1,1,readout_size]);
CS_pattern=squeeze(logical(permute(CS_pattern,[4,2,3,1])));
close all

mask=CS_pattern;

%save image
mask2figure(K_pattern,file_name_img);
save(file_name_mat,'pdf','cs_factor','mask');
end
