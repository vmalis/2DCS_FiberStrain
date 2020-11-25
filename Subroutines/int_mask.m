function [I_mask]=int_mask(magnitude)

% creates mask for velocity images based on mag intensities
%   
% Input:   
%      magnitude data
%
% Output:  
%      mask
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis

I_mask=zeros(size(magnitude));

for i=1:size(magnitude,3)

      a=abs(magnitude(:,:,i));
      thr=mean(a(:))/3.2;
      a(a<thr)=0;
      a(a>=thr)=1;
      se = strel('disk',1); 
      b = imerode(a,se);  
      b = bwareaopen(b,800);
      I_mask(:,:,i)=imfill(b,'holes');
      
end



