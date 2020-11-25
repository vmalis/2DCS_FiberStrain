function [data] = vepc_simple_recon(kspace,CM)

% reconstruction of multidimensional kspace data
%   
% Input:   
%      kspace data
%
% Output:  
%      array of reconstructed velocities and magnitude
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis


[nx,ny,nt,nv,nc]=size(kspace);
Recon_temp=zeros(nx,ny,nt,nv,nc);
Data=zeros(nx,ny,nt,nv);

%% recon from fourier space

for t=1:nt
    for v=1:nv
        for c=1:nc
            %Recon_temp(:,:,t,v,c)=fftshift(fft2(kspace(:,:,t,v,c)),1);
            Recon_temp(:,:,t,v,c)=ifft2c_mri(kspace(:,:,t,v,c));
        end
    end
end

%% extract magnitude and velocity data-------------------------------------

% combine channels
for t=1:nt
    for v=1:nv
        Data(:,:,t,v)=sum(squeeze(Recon_temp(:,:,t,v,:)).*conj(CM),3)./sum(abs(CM).^2,3);
    end
end

% mag & velocity
Data_magnitude=abs(sum(Data,4));
for frame=1:size(Data_magnitude,3)
     Data_magnitude(:,:,frame)=mat2gray(Data_magnitude(:,:,frame));
end

Data_phase = angle(Data);
[v_x,v_y,v_z]=vepc_calc(Data_phase,Data_magnitude);

size(v_x)
size(v_y)
size(v_z)
size(Data_magnitude)

data=cat(4,v_x,v_z,v_y,Data_magnitude);

end 