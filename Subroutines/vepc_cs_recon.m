function [recon_data] = vepc_cs_recon(kspace,mask,CM_w,niteFT,nitePCA,lambda_fft,lambda_pca)


% reconstruction of multidimensional kspace data,
% each phase data seperatly
%   
% Input:   
%      kspace data
%      coil maps
%      number of iterations FT
%      number of iterations PCA
%      lagrange multiplier
%
% Output:  
%      array of reconstructed velocities and magnitude
%
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis

%normalizing data and get cs pattern from cs mask

kspace = kspace/max(max(max(abs(ifft2c(kspace))))) + eps;

cs_pattern=squeeze(mask(:,:,:,1,1));
% k-t sparse sense recon parameters
param.E = Emat_xyt(cs_pattern,CM_w);
param.W = TempFFT(3);
param.TV = TVOP();
param.TVWeight=0;
param.nite = 10;
param.display=1;


tic %start timer

if size(size(kspace),2)>4
  
    recon_data=zeros(size(kspace,1), size(kspace,2), size(kspace,3), size(kspace,4));
    
    for echo=1:size(kspace,4)
        
        kdata=squeeze(kspace(:,:,:,echo,:));
        param.y = kdata;
        
        % linear reconstruction
        recon_dft=param.E'*kdata;
        % ktsparsesense
        recon_data(:,:,:,echo) = ktsparsesense(recon_dft,param,lambda_fft,lambda_pca,niteFT,nitePCA);
         
    end
  
else
   
    recon_dft=param.E'*kspace;
    param.y=kspace;
    recon_data=ktsparsesense(recon_dft,param,lambda_fft,lambda_pca,niteFT,nitePCA);
    
end

toc %stop timer


end


function recon_dft = ktsparsesense(recon_dft,param,lambda_fft,lambda_pca,niteFT,nitePCA)
% k-t SPARSE-SENSE reconstruction

        param.L1Weight=lambda_fft;
        param.L1Weight_pca=0;
        param.nite_fft=niteFT;
        param.nite_pca=0;
        fprintf('\n Reconstruction using temporal FFT');
        fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
        % tempFFT
        for n=1:2
            recon_dft = CSL1NlCg(recon_dft,param);
        end 
        
        
        if nitePCA > 0
            param.L1Weight=0;
            param.L1Weight_pca=lambda_pca;
            param.nite_fft=0;
            param.nite_pca=nitePCA;
            fprintf('\n Reconstruction using temporal PCA');
            fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
            % tempPCA
            for n=1:2
                recon_dft = CSL1NlCg(recon_dft,param);
            end 
        else
        end
end