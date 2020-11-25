%% -----------------------------------------------------------------------   
%                       CS reconstruction for VEPC
%
%  UC San Diego / March 2019 / Vadim Malis
%
%  description: Top level script 1
%
% ------------------------------------------------------------------------
clc
clear all

disp("Starting reconstruction")
disp("=======================================")

%% read pfile
disp("Step 1 of 10: reading pfile to sturcture")
disp("=======================================")
raw=pfile2struct;



for files=1:size(raw,2)


%% may require deletion of frames
% a = [1];
% a=[1,2:2:15];         % cs: 4, 2    vps: 4
% a=[1,2:2:17];         % cs: 3, 2.4  vps: 4
%  a=[1,2:2:31];          % cs: 4       vps: 2  
% % a=[1,2:2:35];         % cs: 3, 2.4  vps: 2
% % %%% remove
% raw(files).kspace(:,:,a,:,:)=[];
% clearvars a


%% read mask or create mask
disp("Step 2 of 10: reading/generating undersampling mask")
disp("=======================================")

% two possible cases: 
%       - read mask from already undersampled data  (get_csmask)
%       - generate mask                             (gen_csmask)
%
% in both cases the mask is saved as a png as well as pattern and mask in
% struct with info

% Get raw data sizes
%[nreadout,nphases,nframes,nechos,ncoils]=size(raw(files).kspace);

% % CASE 1:
cs_factor=3;
[nreadout,nphases, nframes,nechos,ncoils]=size(raw(1).kspace);

[~,mask,~]=gen_csmask(nphases,nreadout,nframes,cs_factor);
clearvars cs_factor

% % CASE 2: or get_csmask()
% [mask]=get_csmask(raw(files).kspace);

% replicate to make same size as the k-space data
mask=repmat(mask,[1,1,1,nechos,ncoils]);

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


%% fully sampled data recon
disp("Step 4 of 10: fully sampled data recon")
disp("=======================================")

data=vepc_simple_recon(raw(files).kspace,CM_w);
suffix=[raw(files).series(1:3),'fs'];

% plot
[Vx,Vy,Vz,M]=vepc_plot(data,raw(files),suffix,1);
n=1;

% fill up struture---------------------
Recon(n).type       =   'original';
Recon(n).L1w_fft    =   NaN;
Recon(n).nite_fft   =   NaN;
Recon(n).L1w_pca    =   NaN;
Recon(n).nite_pca   =   NaN;
Recon(n).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
%   cropped velocities, best for montage
Recon(n).M          =   M;
Recon(n).Vx         =   Vx;
Recon(n).Vy         =   Vy;
Recon(n).Vz         =   Vz;
n=n+1;
% --------------------------------------

%% zero filled CS recon
disp("Step 5 of 10: : zero filled data recon")
disp("=======================================")

data=vepc_simple_recon(raw(files).kspace.*double(mask),CM_w);
suffix=[raw(files).series(1:3),'zf'];

%plot
[Vx,Vy,Vz,M]=vepc_plot(data,raw(files),suffix,1);

%fill up struture---------------------
Recon(n).type       =   'zero filled';
Recon(n).L1w_fft    =   NaN;
Recon(n).nite_fft   =   NaN;
Recon(n).L1w_pca    =   NaN;
Recon(n).nite_pca   =   NaN;
Recon(n).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
  % cropped velocities, best for montage
Recon(n).M          =   M;
Recon(n).Vx         =   Vx;
Recon(n).Vy         =   Vy;
Recon(n).Vz         =   Vz;
n=n+1;
%--------------------------------------


%% CS recon split
disp("Step 6 of 10: : cs recon each phase set")
disp("=======================================")


% nite_fft = [10,     20,     30];
% nite_pca = [0,      10,     20,      30];
% L1_fft   = [0.01,   0.025,  0.05,    0.075,   0.1,    0.15];
% L1_pca   = [0,      0.01,   0.025,   0.05,    0.075,  0.1,      0.15];

% 
nite_fft = 20;
nite_pca = 30;
L1_fft   = 0.075;
L1_pca = 0.075;


temp=combvec(nite_fft,nite_pca,L1_fft,L1_pca);

a=[];

for i = 1:size(temp,2)
    if temp(2,i)==0 && temp(4,i)>0 
        a=[a, i];
    elseif temp(2,i)>0 && temp(4,i)==0
        a=[a, i];
    elseif temp(4,i)>0 && temp(3,i) ~= temp(4,i)
        a=[a, i];
    end
end

a=unique(a);
temp(:,a)=[];

nite_fft = temp(1,:);
nite_pca = temp(2,:);
L1_fft   = temp(3,:);
L1_pca   = temp(4,:);

clearvars temp

wb = waitbar(0,'Not-joint recon: 0% complete');
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
    
%     v_x=cat(3,v_x(:,:,1),v_x(:,:,1)/2,v_x);
%     v_y=cat(3,v_y(:,:,1),v_y(:,:,1)/2,v_y);
%     v_z=cat(3,v_z(:,:,1),v_z(:,:,1)/2,v_z);
%     Data_magnitude=cat(3,Data_magnitude(:,:,1:2),Data_magnitude);
    
    data=cat(4,v_x,v_y,v_z,Data_magnitude);
    
    %plot
    [Vx,Vy,Vz,M]=vepc_plot(data,raw(files),suffix,files);

    %fill up struture---------------------
    Recon(n).series_num =   raw(files).series(1:2);
    name_id                 =   pwd;
    Recon(n).ID         =   name_id(end-1:end);
    Recon(n).location   =   raw(files).header.SeriesData.start_loc;
    Recon(n).MVC        =   raw(files).series(7:8);
    Recon(n).type       =   'not-joint';
    Recon(n).L1w_fft    =   L1_fft(m);
    Recon(n).nite_fft   =   nite_fft(m);
    Recon(n).L1w_pca    =   L1_pca(m);
    Recon(n).nite_pca   =   nite_pca(m);
    Recon(n).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
      % cropped velocities, best for montage
    Recon(n).M          =   M;
    Recon(n).Vx         =   Vx;
    Recon(n).Vy         =   Vy;
    Recon(n).Vz         =   Vz;
    
    
   
    
    n=n+1;
    %--------------------------------------
    frac = m/p;
    msg = sprintf('Not-joint recon: %d%% complete',round(frac*100));
    waitbar(m/p,wb,msg)
    
end
close(wb)


% %% CS recon joint  ref + velocity
% disp("Step 7 of 10: : cs recon joint ref + velocity")
% disp("=======================================")
% 
% %
% nite_fft = [10,     20,     30];
% nite_pca = [0,      10,     20,      30];
% L1_fft   = [0.01,   0.025,  0.05,    0.075,   0.1,    0.15];
% L1_pca   = [0,      0.01,   0.025,   0.05,    0.075,  0.1,      0.15];
% 
% 
% temp=combvec(nite_fft,nite_pca,L1_fft,L1_pca);
% 
% a=[];
% 
% for i = 1:size(temp,2)
%     if temp(2,i)==0 && temp(4,i)>0 
%         a=[a, i];
%     elseif temp(2,i)>0 && temp(4,i)==0
%         a=[a, i];
%     elseif temp(4,i)>0 && temp(3,i) ~= temp(4,i)
%         a=[a, i];
%     end
% end
% 
% a=unique(a);
% temp(:,a)=[];
% 
% nite_fft = temp(1,:);
% nite_pca = temp(2,:);
% L1_fft   = temp(3,:);
% L1_pca   = temp(4,:);
% 
% clearvars temp
% 
% 
% wb = waitbar(0,'Joint recon ref + velocity: 0% complete');
% p=size(nite_fft,2);
% 
% % join data
% temp=raw.kspace.*double(mask);
% data_joint_x = reshape(permute(temp(:,:,:,[1,2],:),[1,2,4,3,5]),[nreadout,nphases,nframes*2,ncoils]);
% data_joint_y = reshape(permute(temp(:,:,:,[1,3],:),[1,2,4,3,5]),[nreadout,nphases,nframes*2,ncoils]);
% data_joint_z = reshape(permute(temp(:,:,:,[1,4],:),[1,2,4,3,5]),[nreadout,nphases,nframes*2,ncoils]);
% mask_joint   = reshape(permute(mask(:,:,:,[1,2],:),[1,2,4,3,5]),[nreadout,nphases,nframes*2,ncoils]);
% 
% clearvars temp
% 
% 
% for m=1:p
%     
%     data_x=vepc_cs_recon(data_joint_x,mask_joint,CM,nite_fft(m),nite_pca(m),L1_fft(m),L1_pca(m));
%     data_y=vepc_cs_recon(data_joint_y,mask_joint,CM,nite_fft(m),nite_pca(m),L1_fft(m),L1_pca(m));
%     data_z=vepc_cs_recon(data_joint_z,mask_joint,CM,nite_fft(m),nite_pca(m),L1_fft(m),L1_pca(m));
%     
%     %unjoin
%     data_x=reshape(data_x,[nreadout,nphases,nframes,2]);
%     data_y=reshape(data_y,[nreadout,nphases,nframes,2]);
%     data_z=reshape(data_z,[nreadout,nphases,nframes,2]);
%     
%     data=cat(4,mean(cat(4,data_x(:,:,:,1),data_y(:,:,:,1),data_z(:,:,:,1)),4),data_x(:,:,:,2),data_y(:,:,:,2),data_z(:,:,:,2));
%     
%    
%     L1_fft_str   = num2str(L1_fft(m));
%     L1_fft_str   = strrep(L1_fft_str,'.','');
%     L1_pca_str   = num2str(L1_pca(m));
%     L1_pca_str   = strrep(L1_pca_str,'.','');
%     
%     suffix=['cssj_' num2str(nite_fft(m)) '_' num2str(nite_pca(m)) '_' L1_fft_str '_' L1_pca_str];
%     
%     
%     % mag & velocity calculation
%     Data_magnitude=abs(sum(data,4));
%     for frame=1:size(Data_magnitude,3)
%         Data_magnitude(:,:,frame)=mat2gray(Data_magnitude(:,:,frame));
%     end
%     Data_phase = angle(data);
%     [v_x,v_y,v_z]=vepc_calc(Data_phase,Data_magnitude);
%     data=cat(4,v_x,v_y,v_z,Data_magnitude);
%     
%     %plot
%     [Vx,Vy,Vz,M]=vepc_plot(data,raw,suffix);
% 
%     %fill up struture---------------------
%     Recon(n).type       =   'incomplete joint';
%     Recon(n).L1w_fft    =   L1_fft(m);
%     Recon(n).nite_fft   =   nite_fft(m);
%     Recon(n).L1w_pca    =   L1_pca(m);
%     Recon(n).nite_pca   =   nite_pca(m);
%     Recon(n).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
%       % cropped velocities, best for montage
%     Recon(n).M          =   M;
%     Recon(n).Vx         =   Vx;
%     Recon(n).Vy         =   Vy;
%     Recon(n).Vz         =   Vz;
%     n=n+1;
%     %--------------------------------------
%     frac = m/p;
%     msg = sprintf('Incomplete joint recon: %d%% complete',round(frac*100));
%     waitbar(m/p,wb,msg)
%     
% end
% close(wb)
% 
% 
% %% CS recon joint
% disp("Step 8 of 10: : cs recon joint all")
% disp("=======================================")
% 
% 
% nite_fft = [10,     20,     30];
% nite_pca = [0,      10,     20,      30];
% L1_fft   = [0.005,  0.01,   0.025,    0.05,    0.075,   0.1];
% L1_pca   = [0,      0.005,   0.01,    0.025,   0.05,    0.075, 0.1];
% 
% 
% temp=combvec(nite_fft,nite_pca,L1_fft,L1_pca);
% 
% a=[];
% 
% for i = 1:size(temp,2)
%     if temp(2,i)==0 && temp(4,i)>0 
%         a=[a, i];
%     elseif temp(2,i)>0 && temp(4,i)==0
%         a=[a, i];
%     elseif temp(4,i)>0 && temp(3,i) ~= temp(4,i)
%         a=[a, i];
%     end
% end
% 
% a=unique(a);
% temp(:,a)=[];
% 
% nite_fft = temp(1,:);
% nite_pca = temp(2,:);
% L1_fft   = temp(3,:);
% L1_pca   = temp(4,:);
% 
% clearvars temp
% 
% 
% wb = waitbar(0,'Joint recon: 0% complete');
% p=size(nite_fft,2);
% 
% % join data
% data_joint=reshape(permute(raw.kspace.*double(mask),[1,2,4,3,5]),[nreadout,nphases,nframes*nechos,ncoils]);
% mask_joint=reshape(permute(mask,[1,2,4,3,5]),[nreadout,nphases,nframes*nechos,ncoils]);
% 
% for m=29:p
%     
%     data=vepc_cs_recon(data_joint,mask_joint,CM,nite_fft(m),nite_pca(m),L1_fft(m),L1_pca(m));
%     
%     %unjoin
%     data=reshape(data,[nreadout,nphases,nframes,nechos]);
%     
%     L1_fft_str   = num2str(L1_fft(m));
%     L1_fft_str   = strrep(L1_fft_str,'.','');
%     L1_pca_str   = num2str(L1_pca(m));
%     L1_pca_str   = strrep(L1_pca_str,'.','');
%     
%     suffix=['csj_' num2str(nite_fft(m)) '_' num2str(nite_pca(m)) '_' L1_fft_str '_' L1_pca_str];
%     
%     
%     % mag & velocity calculation
%     Data_magnitude=abs(sum(data,4));
%     for frame=1:size(Data_magnitude,3)
%         Data_magnitude(:,:,frame)=mat2gray(Data_magnitude(:,:,frame));
%     end
%     Data_phase = angle(data);
%     [v_x,v_y,v_z]=vepc_calc(Data_phase,Data_magnitude);
%     data=cat(4,v_x,v_y,v_z,Data_magnitude);
%     
%     %plot
%     [Vx,Vy,Vz,M]=vepc_plot(data,raw,suffix);
% 
%     %fill up struture---------------------
%     Recon(n).type       =   'joint';
%     Recon(n).L1w_fft    =   L1_fft(m);
%     Recon(n).nite_fft   =   nite_fft(m);
%     Recon(n).L1w_pca    =   L1_pca(m);
%     Recon(n).nite_pca   =   nite_pca(m);
%     Recon(n).data       =   data;   % x,y,frame, [mag, vx, vy,vz]
%       % cropped velocities, best for montage
%     Recon(n).M          =   M;
%     Recon(n).Vx         =   Vx;
%     Recon(n).Vy         =   Vy;
%     Recon(n).Vz         =   Vz;
%     n=n+1;
%     %--------------------------------------
%     frac = m/p;
%     msg = sprintf('Joint recon: %d%% complete',round(frac*100));
%     waitbar(m/p,wb,msg)
%     
% end
% close(wb)


% move images to appropriate folders
disp("Step 8 of 9: : cs recon")
disp("=======================================")

mkdir('mag')
file=dir('mag*.png');
wb = waitbar(0,'moving magnitude images: 0% complete');
p=size(file,1);

for m=1:p
    movefile(file(m).name, 'mag')
    frac = m/p;
    msg = sprintf('moving magnitude images: %d%% complete',round(frac*100));
    waitbar(m/p,wb,msg)
end
close(wb)

mkdir('vel')
file=dir('vel*.png');
wb = waitbar(0,'moving velocity images: 0% complete');
p=size(file,1);

for m=1:p
    movefile(file(m).name, 'mag')
    frac = m/p;
    msg = sprintf('moving velocity images: %d%% complete',round(frac*100));
    waitbar(m/p,wb,msg)
end
close(wb)





%% result comparsion
disp("Step 9 of 9: : cs recon")
disp("=======================================")
fprintf("\n Error calculations \n");

%original data smoothing

Vx_ftrd=zeros(size(Recon(1).Vx));
Vy_ftrd=Vx_ftrd;
Vz_ftrd=Vx_ftrd;


for frame=1:nframes
    Vx_ftrd(:,:,frame) = anisodiff2D(Recon(1).Vx(:,:,frame),Recon(1).M(:,:,frame),15,1/7,6,1);
    Vy_ftrd(:,:,frame) = anisodiff2D(Recon(1).Vy(:,:,frame),Recon(1).M(:,:,frame),15,1/7,6,1);
    Vz_ftrd(:,:,frame) = anisodiff2D(Recon(1).Vz(:,:,frame),Recon(1).M(:,:,frame),15,1/7,6,1);
end


for i=1:size(Recon,2)
    
    i
    
    %processed data smoothing
    Vx_ftrd_cs=zeros(size(Recon(1).Vx));
    Vy_ftrd_cs=Vx_ftrd;
    Vz_ftrd_cs=Vx_ftrd;
    
    for frame=1:nframes
        Vx_ftrd_cs(:,:,frame) = anisodiff2D(Recon(i).Vx(:,:,frame),Recon(i).M(:,:,frame),15,1/7,6,1);
        Vy_ftrd_cs(:,:,frame) = anisodiff2D(Recon(i).Vy(:,:,frame),Recon(i).M(:,:,frame),15,1/7,6,1);
        Vz_ftrd_cs(:,:,frame) = anisodiff2D(Recon(i).Vz(:,:,frame),Recon(i).M(:,:,frame),15,1/7,6,1);
    end
    
    
    Recon(i).vel_error_x = abs((Vx_ftrd-Vx_ftrd_cs))./abs(Vx_ftrd)*30;
    Recon(i).vel_error_y = abs((Vy_ftrd-Vy_ftrd_cs))./abs(Vy_ftrd)*30;
    Recon(i).vel_error_z = abs((Vz_ftrd-Vz_ftrd_cs))./abs(Vz_ftrd)*30;
    
    Recon(i).vel_error_x(Recon(1).M<0.08)=0;
    Recon(i).vel_error_y(Recon(1).M<0.08)=0;
    Recon(i).vel_error_z(Recon(1).M<0.08)=0;
    
    
    if strcmp(Recon(i).type,'original')
        type='o';
        elseif strcmp(Recon(i).type,'zero filled')
        type='zf';
        elseif strcmp(Recon(i).type,'not-joint')
        type='nj';
        elseif strcmp(Recon(i).type,'incomplete joint')
        type='sj';
        else strcmp(Recon(i).type,'joint')
        type='j';
    end
        
    strrep(num2str(Recon(i).nite_fft),'.','');
    
    name=['dv_',type,'_',num2str(Recon(i).nite_fft),'_',...
            strrep(num2str(Recon(i).L1w_fft),'.',''),'_',...
            num2str(Recon(i).nite_pca),'_',...
            strrep(num2str(Recon(i).L1w_pca),'.','')];
    %plot
    vepc_error_plot(Recon(i).vel_error_x,Recon(i).vel_error_y,Recon(i).vel_error_z,name)
    
    
end

save Recon.mat Recon -v7.3



end




