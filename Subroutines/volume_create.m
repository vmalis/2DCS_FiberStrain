%%% merge volumes, do calcs plots colormaps and stats

clc
clear all

%% volume merge
load('force_data.mat')
load('cs_processed.mat')
cd DTI
load('DTI.mat')
cd ..

% important parameters
MVC=unique({Recon.MVC},'stable');
num_iter=10; kappa=6; option=1; delta_t=1/7;
pixelspacing = [1.1719,1.1719,5.0000]';


for volume=1:size(unique({Recon.MVC}),2)

    Data(volume).ID = unique([Recon.ID],'stable');
    Data(volume).MVC= MVC(volume);
    
    same_mvc_series = find(contains({Recon.MVC},MVC(volume)));
    
    [~,~,dti_locations] = intersect(round([Recon(same_mvc_series).location]),round(DTI_data.location'),'stable');
    
    Data(volume).vx = permute(cat(4,Recon(same_mvc_series(1)).Vx,Recon(same_mvc_series(2)).Vx,Recon(same_mvc_series(3)).Vx),[1,2,4,3]);
    Data(volume).vy = permute(cat(4,Recon(same_mvc_series(1)).Vy,Recon(same_mvc_series(2)).Vy,Recon(same_mvc_series(3)).Vy),[1,2,4,3]);
    Data(volume).vz = permute(cat(4,Recon(same_mvc_series(1)).Vz,Recon(same_mvc_series(2)).Vz,Recon(same_mvc_series(3)).Vz),[1,2,4,3]);
    
    mask=zeros(size(Data(volume).vx));
    mask(Data(volume).vx~=0)=1;
    Data(volume).mask = mask;
    Data(volume).m  = mask.*permute(cat(4,Recon(same_mvc_series(1)).M,Recon(same_mvc_series(2)).M,Recon(same_mvc_series(3)).M),[1,2,4,3]);
    
    Data(volume).dti_vec  = DTI_data.DTI_eigenvector(:,:,dti_locations,:,:);
    Data(volume).dti_val  = DTI_data.DTI_eigenvalue(:,:,dti_locations,:);
    
    [~,~,force_series] = intersect({Recon(same_mvc_series).series_num},{Force.series_num},'stable');
    
    Data(volume).MVC_true = mean([Force(force_series).pcent]);
    Data(volume).force    = mean(reshape([Force(force_series).mean],[size(Force(force_series(1)).mean,2),size(force_series,1)]),2);
    
    frames=size(Data(volume).vx,4);
    
    
    h=waitbar(0,'Filtering');
    
    for k=1:frames

         Data(volume).vx_sm(:,:,:,k) = anisodiff3D(Data(volume).vx(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);
         Data(volume).vy_sm(:,:,:,k) = anisodiff3D(Data(volume).vy(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);
         Data(volume).vz_sm(:,:,:,k) = anisodiff3D(Data(volume).vz(:,:,:,k),...
                              num_iter,delta_t,kappa,option,pixelspacing);

         waitbar(k/frames,h)
         
    end
    
    close(h)
end

clearvars -except Data frames pixelspacing
close all

%% ROI and calcs


a=figure;
imshow(mat2gray(Data(series).m(:,:,2,1)),'Initialmagnification',400)
h_all = drawfreehand(gca);
mask_all = roiWait(h_all);
close(a)

a=figure;
imshow(mat2gray(Data(series).m(:,:,2,1)),'Initialmagnification',400);
h_gm  = drawrectangle(gca,'Position',[50 50 10 30]);
mask_gm = roiWait(h_gm);
close(a)

a=figure;
imshow(mat2gray(Data(series).m(:,:,2,1)),'Initialmagnification',400);
h_sol  = drawrectangle(gca,'Position',[50 50 10 30]);
mask_sol = roiWait(h_sol);
close(a)


for series=1:size(Data,2)

%predifined from acqusition
trigger=ones(frames-1,1)*0.136;


% create coordinates

mask=Data(series).mask(:,:,:,:);
mask(mask==0)=NaN;
[x,y,z]=meshgrid(1:256,1:256,1:3);

x(isnan(mask(:,:,:,1)))=NaN;
y(isnan(mask(:,:,:,1)))=NaN;
z(isnan(mask(:,:,:,1)))=NaN;

x_nan=x(:);
y_nan=y(:);
z_nan=z(:);

x_nan(isnan(x_nan))=[];
y_nan(isnan(y_nan))=[];
z_nan(isnan(z_nan))=[];



% tracking
%
% starting from here the velocities are in image plane convention

VX_sm=Data(series).vx_sm.*mask;
VY_sm=Data(series).vz_sm.*mask;
VZ_sm=Data(series).vy_sm.*mask;

tic
[xs,ys,zs,vr,vx,vy,vz] = track3d(x_nan,y_nan,z_nan,VX_sm,VY_sm,VZ_sm,trigger,pixelspacing);  
toc                 
                                
% obtain tracked and padd with zeros needed cause we dont track nans
% big time saver

XS=zeros(size(Data(series).vx));
YS=XS;
ZS=XS;
VR=XS;
VX=XS;
VY=XS;
VZ=XS;


for t=1:frames
    for index=1:size(x_nan)
        XS(y_nan(index),x_nan(index),z_nan(index),t)=xs(index,t);
        YS(y_nan(index),x_nan(index),z_nan(index),t)=ys(index,t);
        ZS(y_nan(index),x_nan(index),z_nan(index),t)=zs(index,t);
        VR(y_nan(index),x_nan(index),z_nan(index),t)=vr(index,t);
        VX(y_nan(index),x_nan(index),z_nan(index),t)=vx(index,t);
        VY(y_nan(index),x_nan(index),z_nan(index),t)=vy(index,t);
        VZ(y_nan(index),x_nan(index),z_nan(index),t)=vz(index,t);
    end
end

% time dependent mask
mask_dynamic      = dynamic_mask(XS,YS,ZS,squeeze(mask(:,:,1,1)));
mask_all_dynamic  = dynamic_mask(XS,YS,ZS,mask_all);
mask_sol_dynamic  = dynamic_mask(XS,YS,ZS,mask_sol);
mask_gm_dynamic   = dynamic_mask(XS,YS,ZS,mask_gm);
                                
                              
% tensor calcs                                
Data(series).strain     =   strain_tensor_calcs(VX,VY,VZ,XS,YS,ZS,pixelspacing,trigger);
[Data(series).strain.ED, Data(series).strain.LD, Data(series).strain.SRD]  =   fa_strain_calcs(Data(series).strain,Data(series).dti_vec);

% masking
Data(series).strain_SOL = strain_ROI(Data(series).strain,double(mask_sol_dynamic));
Data(series).strain_GM  = strain_ROI(Data(series).strain,double(mask_gm_dynamic));

% peak
Data(series).strain_SOL_peak = peak_ROI(Data(series).strain_SOL);
Data(series).strain_GM_peak  = peak_ROI(Data(series).strain_GM);


%colormaps
foldername=sprintf(Data(series).MVC);
mkdir(foldername)
cd(foldername)
strain_colormaps(Data(series).strain)
cd ..



end


mkdir('plots')
cd('plots')

% GM

strain_mvc_plots(Data.strain_GM,'GM')

% SOL
strain_mvc_plots(Data.strain_SOL,'SOL')

save('Results.mat',Data)