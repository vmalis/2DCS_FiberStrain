%% simple rigid body registration between two images

function [moving_registered]=mri_register(fixed,moving)

moving_registered=zeros(size(fixed));

    for slice=1:size(fixed,3)

        [optimizer,metric] = imregconfig('multimodal');
        optimizer.MaximumIterations = 1000;
        correctedI = imregister(squeeze(moving(:,:,slice)),squeeze(fixed(:,:,slice)),'affine',optimizer,metric);
        moving_registered(:,:,slice)=correctedI;
    end



end
