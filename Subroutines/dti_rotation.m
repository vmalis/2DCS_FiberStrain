%=========================================
%   Step4 for 3d Strain Rate Analysis
%=========================================
%
% Final scipt to extract peak values and rotate SR to DTI
%_____________________________________________________
% written by Vadim Malis
% 10/17 at UCSD RIL

function [SR_D]=dti_rotation(rotation_matrix,tensor_to_rotate)

SR_D=zeros(size(tensor_to_rotate));


ii=size(tensor_to_rotate,1);
jj=size(tensor_to_rotate,2);
kk=size(tensor_to_rotate,3);
tt=size(tensor_to_rotate,4);

    for i=1:ii
        for j=1:jj
            for k=1:kk
                for t=1:tt
                                        
                    R=squeeze(rotation_matrix(i,j,k,:,:));
                    Q=[R(:,3),R(:,2),R(:,1)];
                    SR_D(i,j,k,t,:,:)=Q*...
                                squeeze(tensor_to_rotate(i,j,k,t,:,:))*...
                                Q';

                     end
                    
                end
            end
        end
    end
    
