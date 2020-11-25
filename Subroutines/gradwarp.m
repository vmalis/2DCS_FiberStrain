%% Function to perform geometric corrections using GE gradwarp
%
% Important:  original pfile must be loaded and set to active

function [data_unwarped] = gradwarp(data_warped,corners)

    for i=1:size(data_warped,3)

        temp=imresize(data_warped(:,:,i),[max(size(data_warped)),max(size(data_warped))]);
        data_unwarped(:,:,i) = GERecon('Gradwarp', temp, corners);
        
    end
     
end