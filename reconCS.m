function Recon = reconCS(raw,type,N)
%% -----------------------------------------------------------------------
%  UC San Diego / November 2020 / Vadim Malis
%
%  Subroutine to read raw data and perform recon
%
%   input:  raw data structure -> from pfile2struct
%           type:       'full' or 'cs*vps*' 
%           Number:     integer
%-------------------------------------------------------------------------


switch type
    
    case 'full'
        
        a = [1];
        disp('full sampling recon')
        
    case 'cs4vps4'
        
        disp('CS-factor = 4, VPS = 4')
        a=[1,2:2:15];
        
    case 'cs2vps4'
        
        disp('CS-factor = 4, VPS = 4')
        a=[1,2:2:15];
        
    case 'cs3vps4'
        
        disp('CS-factor = 3, VPS = 4')
        a=[1,2:2:17];   
        
    case 'cs2.4vps4'
        
        disp('CS-factor = 2.4, VPS = 4')
        a=[1,2:2:17];   
        
    case 'cs4vps2'
        
        a=[1,2:2:31];
        disp('CS-factor = 4, VPS = 2')
        
    case 'cs3vps2'
        
        a=[1,2:2:35];
        disp('CS-factor = 3, VPS = 2')
        
    case 'cs2.4vps2'
        
        a=[1,2:2:35]; 
        disp('CS-factor = 2.4, VPS = 2')  
        
    otherwise
        disp('Unknown case')
        return
    
end


%% cut not sampled frames
raw.kspace(:,:,a,:,:)=[];
clearvars a

%% read undersampling mask
disp("Reading/generating undersampling mask")
disp("=======================================")

% Get raw data sizes
[nreadout,nphases,nframes,nechos,ncoils]=size(raw.kspace);

% Get_csmask()
[mask]=get_csmask(raw.kspace);

% Replicate to make same size as the k-space data
mask=repmat(mask,[1,1,1,nechos,ncoils]);



%% create coil maps
disp("Obtaining coilmap")
disp("=======================================")

%-----adaptive array recon if multichannel-----
if raw.channels>1
    temp=squeeze(mean(mean(raw.kspace,3),4));
    Temp=zeros(size(temp));
    
    for coil=1:size(Temp,3)
        %Temp(:,:,coil)=fftshift(fft2(temp(:,:,coil)),1);    
        Temp(:,:,coil)=ifft2c_mri(temp(:,:,coil));    
        
    end
    
    [I,CM]=adapt_array_cmap(Temp,eye(ncoils),1);
    %CM=flip(CM,1);
    coilmap2figure(CM,raw,N)
else
    CM=ones(raw.xres,raw.yres);
end

%--------coil weights from the prescan--------- 
info=GERecon('Pfile.Header');
coil_weights=info.PrescanHeader.rec_std(1:ncoils);
CM_w=zeros(size(CM));

for i=1:ncoils
    CM_w(:,:,i)=CM(:,:,i)*(coil_weights(i));
end

clearvars coil temp Temp I


switch type
    
    case 'full'
        
        disp("Full k-space recon")
        disp("=======================================")
        
        data=vepc_simple_recon(raw.kspace,CM_w);
        data=flip(data,1);
        suffix=[raw.series(1:3),'fs'];

        % plot
        
        [Vx,Vy,Vz,M]=vepc_plot(data,raw,suffix,N);

        % fill up struture---------------------
        Recon.type       =   'original';
        Recon.L1w_fft    =   NaN;
        Recon.nite_fft   =   NaN;
        Recon.L1w_pca    =   NaN;
        Recon.nite_pca   =   NaN;
        Recon.data       =   data;   % x,y,frame, [mag, vx, vy,vz]
        %   cropped velocities, best for montage
        Recon.M          =   M;
        Recon.Vx         =   Vx;
        Recon.Vy         =   Vy;
        Recon.Vz         =   Vz;
        Recon.ID         =   pwd;

    otherwise
        
        disp("CS recon")
        disp("=======================================")

        nite_fft = 20;
        nite_pca = 30;
        L1_fft   = 0.025;
        L1_pca = 0.025;

        clearvars temp

        p=size(nite_fft,2);

        for m=1:p
    
            data=vepc_cs_recon(raw.kspace.*double(mask),mask,CM,nite_fft(m),nite_pca(m),L1_fft(m),L1_pca(m));
            data=flip(data,1);
            L1_fft_str   = num2str(L1_fft(m));
            L1_fft_str   = strrep(L1_fft_str,'.','');
            L1_pca_str   = num2str(L1_pca(m));
            L1_pca_str   = strrep(L1_pca_str,'.','');
    
            suffix=[raw.series(1:3), 'cs_' num2str(nite_fft(m)) '_' num2str(nite_pca(m)) '_' L1_fft_str '_' L1_pca_str];
    
            % mag & velocity calculation
            Data_magnitude=abs(sum(data,4));
            for frame=1:size(Data_magnitude,3)
                Data_magnitude(:,:,frame)=mat2gray(Data_magnitude(:,:,frame));
            end
            Data_phase = angle(data);
            
            [v_x,v_y,v_z]=vepc_calc(Data_phase,Data_magnitude);
            v_x=cat(3, v_x(:,:,1)/8, v_x(:,:,1)/2, v_x);
            v_y=cat(3, v_y(:,:,1)/8, v_y(:,:,1)/2, v_y);
            v_z=cat(3, v_z(:,:,1)/8, v_z(:,:,1)/2, v_z);
            
            Data_magnitude=cat(3,Data_magnitude(:,:,1:2),Data_magnitude);

            data=cat(4,v_x,v_y,v_z,Data_magnitude);
    
            %plot
            [Vx,Vy,Vz,M]=vepc_plot(data,raw,suffix,N);

            %fill up struture---------------------
            Recon.series_num =   raw.series(1:2);
            Recon.ID         =   pwd;
            Recon.location   =   raw.header.SeriesData.start_loc;
            Recon.type       =   'not-joint';
            Recon.L1w_fft    =   L1_fft(m);
            Recon.nite_fft   =   nite_fft(m);
            Recon.L1w_pca    =   L1_pca(m);
            Recon.nite_pca   =   nite_pca(m);
            Recon.data       =   data;   % x,y,frame, [mag, vx, vy,vz]
            % cropped velocities, best for montage
            Recon.M          =   M;
            Recon.Vx         =   Vx;
            Recon.Vy         =   Vy;
            Recon.Vz         =   Vz;

        end   

    
        
end
