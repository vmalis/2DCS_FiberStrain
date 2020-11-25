function [E] = tensor_calcs(dx,dy,dz,pixel_spacing)
%   -----------------------------------------------------------------------
%   usage: [xx,xy,xz,yy,yz,zz] = tensor_calcs(dx,dy,dz,res)
%   -----------------------------------------------------------------------
%   input arguments:
%   ----------------
%   dx: x displacement (pxls)       |
%   dy: y displacement (pxls)       | -> all as a 3d array
%   dz: z displacement (slices)     |
%          
%   res: 3 component vector to bring dx,dy,dz to same units
%   -----------------------------------------------------------------------
%   returns:
%   --------
%
%     9 components of a symmetric Lagrangian Strain tensor 
%
%               | Exx   Exy  Exz |
%      E    =   | Eyx   Eyy  Eyz |
%               | Ezx   Ezy  Ezz |
%
%   -----------------------------------------------------------------------
%   written by Vadim Malis
%   11/14 at UCSD RIL
%   -----------------------------------------------------------------------

%identity matrix
I=eye(3);


%deformation matrix
F=zeros(size(dx,1),size(dx,2),size(dx,3),3,3);

[F(:,:,:,1,1),F(:,:,:,1,2),F(:,:,:,1,3)]=gradient(dx,pixel_spacing(1),...
                                                     pixel_spacing(2),...
                                                     pixel_spacing(3));
[F(:,:,:,2,1),F(:,:,:,2,2),F(:,:,:,2,3)]=gradient(dy,pixel_spacing(1),...
                                                     pixel_spacing(2),...
                                                     pixel_spacing(3));
[F(:,:,:,3,1),F(:,:,:,3,2),F(:,:,:,3,3)]=gradient(dz,pixel_spacing(1),...
                                                     pixel_spacing(2),...
                                                     pixel_spacing(3));

%Strain tensor
E=zeros(size(dx,1),size(dx,2),size(dx,3),3,3);

for i=1:size(dx,1)
    for j=1:size(dx,2)
        for k=1:size(dx,3)
            
        f=squeeze(F(i,j,k,:,:));
        E(i,j,k,:,:)=0.5*(f+f');
 
        end
    end
end
end

