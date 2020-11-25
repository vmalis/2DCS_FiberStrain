function im_mean = im_mean(im)
%-------------------------------------------------------------------------%
% im_mean = im_mean(im)
% im is a three-dimensional array
% im_mean returns the two-dimensional mean image
% David Shin
% 2008/02/28
%-------------------------------------------------------------------------%
num_im = size(im,3);

r = size(im,1);
c = size(im,2);

im_sum = zeros(r,c);

for i = 1:num_im
    im_sum = im_sum + im(:,:,i);
end

im_mean = im_sum / num_im;
    

