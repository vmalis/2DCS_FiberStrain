%% flatten the 3d array to 2d array
function [flat_array2d]=array_fltn3d(array)
    
    [size_1,size_2,size_3]=size(array);
    flat_array2d=reshape(array,[size_1*size_2,size_3]);
    
end