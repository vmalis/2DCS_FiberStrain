%% -----------------------------------------------------------------------%   
%                         CS reconstruction for VEPC                      %
%                                                                         %
%                            cs/vps combinations                          %
%                                                                         %
%                   UC San Diego / June 2019 / Vadim Malis                %
%                                                                         %
% ------------------------------------------------------------------------%


clc
clear all

disp("Starting reconstruction")
disp("=======================================")

%% read pfile
disp("Step 1 of 10: reading pfile to sturcture")
disp("=======================================")
raw=pfile2struct;


raw(1).a=[1,20:22];
raw(2).a=[1,2:2:15];
raw(3).a=[1,2:2:15];
raw(4).a=[1,2:2:31];
raw(5).a=[1,2:2:17];
raw(6).a=[1,2:2:17];
raw(7).a=[1,2:2:37];
raw(8).a=[1,2:2:37];

if size(raw,2)>8
    raw(9:end)=[];
end

for i=1:8
    
    raw(i).kspace(:,:,raw(i).a,:,:)=[];
    raw(i).slices=size(raw(i).kspace,3);
    size(raw(i).kspace)
end




%% 
for files=1:size(raw,2)
    
    
    
    % read mask or create mask
    disp("processing data set")
    files
    disp("=======================================")


    % Get raw data sizes
    [nreadout,nphases,nframes,nechos,ncoils]=size(raw(files).kspace);

    
    if files == 1
        mask=ones(nreadout,nphases,nframes,nechos,ncoils);
    else   
        [mask]=get_csmask(raw(files).kspace);
        % replicate to make same size as the k-space data
        mask=repmat(mask,[1,1,1,nechos,ncoils]);
    end
    
    
%% create coil maps
disp("Step 3 of 10: obtaining coilmap")
disp("=======================================")

%-----adaptive array recon if multichannel-----
if raw(files).channels>1
    temp=squeeze(mean(mean(raw(files).kspace,3),4));
    Temp=zeros(size(temp));
    
    for coil=1:size(Temp,3)
        %Temp(:,:,coil)=fftshift(fft2(temp(:,:,coil)),1);    
        Temp(:,:,coil)=ifft2c_mri(temp(:,:,coil));    
        
    end
    
    [I,CM]=adapt_array_cmap(Temp,eye(ncoils),1);
    coilmap2figure(CM,raw(files),files)
    else
    CM=ones(raw(files).xres,raw(files).yres);
end

%--------coil weights from the prescan--------- 
info=GERecon('Pfile.Header');
coil_weights=info.PrescanHeader.rec_std(1:ncoils);
CM_w=zeros(size(CM));

for i=1:ncoils
    CM_w(:,:,i)=CM(:,:,i)*(coil_weights(i));
end

clearvars coil temp Temp I




if files == 1
    
%% fully sampled data recon


data=vepc_simple_recon(raw(files).kspace,CM_w);
suffix=[raw(files).series(1:3),'fs'];

% plot
[Vx,Vy,Vz,M]=vepc_plot(data,raw(files),suffix,files);


% fill up struture---------------------
Recon(files).type       =   'original';
Recon(files).L1w_fft    =   NaN;
Recon(files).nite_fft   =   NaN;
Recon(files).L1w_pca    =   NaN;
Recon(files).nite_pca   =   NaN;
Recon(files).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
%   cropped velocities, best for montage
Recon(files).M          =   M;
Recon(files).Vx         =   Vx;
Recon(files).Vy         =   Vy;
Recon(files).Vz         =   Vz;



else

nite_fft = 20;
nite_pca = 30;
L1_fft   = 0.025;
L1_pca = 0.025;


clearvars temp

%wb = waitbar(0,'Not-joint recon: 0% complete');
p=size(nite_fft,2);

for m=1:p
    
    data=vepc_cs_recon(raw(files).kspace.*double(mask),mask,CM,nite_fft(m),nite_pca(m),L1_fft(m),L1_pca(m));
    
    L1_fft_str   = num2str(L1_fft(m));
    L1_fft_str   = strrep(L1_fft_str,'.','');
    L1_pca_str   = num2str(L1_pca(m));
    L1_pca_str   = strrep(L1_pca_str,'.','');
    
    suffix=[raw(files).series(1:3), 'csnj_' num2str(nite_fft(m)) '_' num2str(nite_pca(m)) '_' L1_fft_str '_' L1_pca_str];
    
    % mag & velocity calculation
    Data_magnitude=abs(sum(data,4));
    for frame=1:size(Data_magnitude,3)
        Data_magnitude(:,:,frame)=mat2gray(Data_magnitude(:,:,frame));
    end
    Data_phase = angle(data);
    [v_x,v_y,v_z]=vepc_calc(Data_phase,Data_magnitude);
    
    v_x=1.05*cat(3,v_x(:,:,1)/6,v_x(:,:,1)/2,v_x);
    v_y=1.05*cat(3,v_y(:,:,1)/6,v_y(:,:,1)/2,v_y);
    v_z=1.05*cat(3,v_z(:,:,1)/6,v_z(:,:,1)/2,v_z);
    Data_magnitude=cat(3,Data_magnitude(:,:,1:2),Data_magnitude);
    
    %%% permute velocities so xy in plane and z through
    data=cat(4,v_x,v_z,v_y,Data_magnitude);
    
    %plot
    [Vx,Vy,Vz,M]=vepc_plot(data,raw(files),suffix,files);

    %fill up struture---------------------
    Recon(files).series_num =   raw(files).series(1:2);
    name_id                 =   pwd;
    Recon(files).ID         =   name_id(end-1:end);
    Recon(files).location   =   raw(files).header.SeriesData.start_loc;
    Recon(files).MVC        =   raw(files).series(7:8);
    Recon(files).type       =   'not-joint';
    Recon(files).L1w_fft    =   L1_fft(m);
    Recon(files).nite_fft   =   nite_fft(m);
    Recon(files).L1w_pca    =   L1_pca(m);
    Recon(files).nite_pca   =   nite_pca(m);
    Recon(files).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
      % cropped velocities, best for montage
    Recon(files).M          =   M;
    Recon(files).Vx         =   Vx;
    Recon(files).Vy         =   Vy;
    Recon(files).Vz         =   Vz;
    
    
   
    
end  
    
    
    
    
end 
 
end

Recon(1).cs     =   1;
Recon(1).vps    =   4;

Recon(2).cs     =   2;
Recon(2).vps    =   4;

Recon(3).cs     =   4;
Recon(3).vps    =   4;

Recon(4).cs     =   4;
Recon(4).vps    =   2;

Recon(5).cs     =   3;
Recon(5).vps    =   4;

Recon(6).cs     =   2.4;
Recon(6).vps    =   4;

Recon(7).cs     =   3;
Recon(7).vps    =   2;

Recon(8).cs     =   2.4;
Recon(8).vps    =   2;

save cs_processed.mat Recon
