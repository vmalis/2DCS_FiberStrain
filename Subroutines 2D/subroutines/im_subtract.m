function im_sub = im_subtract(im_ref, im)
%-------------------------------------------------------------------------%
% im_sub = im_subtract(im_ref, im)
% im_ref is a three-dimensional array of images
% im is an image to be subracted
% im_sub returns a three-dimensional array of the same size as im_ref
% David Shin
% 2008/02/28
%-------------------------------------------------------------------------%
num_im = size(im_ref,3);
r = size(im_ref,1);
c = size(im_ref,2);

im_sub = zeros(r,c,num_im);

for i = 1:num_im
    im_sub(:,:,i) = im_ref(:,:,i) - im;
end

