function [mask]=get_csmask(kspace)

% Generates undersampling pattern
%   
% Input:   
%      kspace
%
% Output:  
%      mask
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis

[nreadout,~,~,~,~]=size(kspace);

%% mask 
temp=squeeze(sum(sum(sum(abs(kspace),1),4),5));
temp(temp>0)=1;

%% montage and save CS pattern
k_pattern=temp';
cs_factor=size(find(k_pattern),1)/size(k_pattern(:),1);


file_name_img = strcat('CS','_x',num2str(1/cs_factor),'.eps');
file_name_mat = strcat('CS','_x',num2str(1/cs_factor),'.mat');

mask=permute(repmat(k_pattern',[1,1,nreadout]),[3,1,2]);

mask2figure(k_pattern,file_name_img);
save(file_name_mat,'cs_factor','mask');


end