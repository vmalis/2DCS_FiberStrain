function [v_x,v_y,v_z,m]=vepc_plot(data,raw,name_suffix,n)

% velocity calculation from phase data
%   
% Input:   
%      data - [vx,vy,vz,magnitude]
%
% Output:  
%      montage of vepc
% ---------------------------------------
% UC San Diego / March 2019 / Vadim Malis


vx=squeeze(data(:,:,:,1));
vy=squeeze(data(:,:,:,2));
vz=squeeze(data(:,:,:,3));
mag=squeeze(data(:,:,:,4));


% %% mask with intensity mask
[I_mask]=int_mask(mag);

vx=vx.*I_mask;
vy=vy.*I_mask;
vz=vz.*I_mask;

%% correct for gradient nonlinearity
file_list=rdir('*.7');
p=GERecon('Pfile.Load',file_list(n).name);
GERecon('Pfile.SetActive', p)

vx_corrected=zeros(size(vx,1),size(vx,1),size(vx,3));
vy_corrected=vx_corrected;
vz_corrected=vx_corrected;
mag_corrected=vx_corrected;


%% move histo

max_vx=max(vx(:));
max_vy=max(vy(:));
max_vz=max(vz(:));

min_vx=min(vx(:));
min_vy=min(vy(:));
min_vz=min(vz(:));

    for frame=1:size(vx,3)
        
        % correct for gradient nonlinearity
        temp=gradwarp(vx(:,:,frame),raw.corners(1));
        vx_corrected(:,:,frame) = temp;
        
        temp=gradwarp(vy(:,:,frame),raw.corners(1));
        vy_corrected(:,:,frame) = temp;
        
        temp=gradwarp(vz(:,:,frame),raw.corners(1));
        vz_corrected(:,:,frame) = temp;
        
        temp=gradwarp(mag(:,:,frame),raw.corners(1));
        mag_corrected(:,:,frame) = temp;
        
    end

%% crop

max_vx_c=max(vx_corrected(:));
max_vy_c=max(vy_corrected(:));
max_vz_c=max(vz_corrected(:));

min_vx_c=min(vx_corrected(:));
min_vy_c=min(vy_corrected(:));
min_vz_c=min(vz_corrected(:));

vx_corrected(vx_corrected>0)=vx_corrected(vx_corrected>0)/max_vx_c*max_vx;
vy_corrected(vy_corrected>0)=vy_corrected(vy_corrected>0)/max_vy_c*max_vy;
vz_corrected(vz_corrected>0)=vz_corrected(vz_corrected>0)/max_vz_c*max_vz;
vx_corrected(vx_corrected<0)=vx_corrected(vx_corrected<0)/min_vx_c*min_vx;
vy_corrected(vy_corrected<0)=vy_corrected(vy_corrected<0)/min_vy_c*min_vy;
vz_corrected(vz_corrected<0)=vz_corrected(vz_corrected<0)/min_vz_c*min_vz;



    t=find(squeeze(abs(mag_corrected(end,:,1))));
    if isempty(t)
        x1=1;
        x2=size(mag_corrected,2);
    else
        x1=t(1);
        x2=t(end);
    end
    
    VX = flip(vx_corrected(:,x1:x2,:),2);
    v_x= flip(vx_corrected(:,:,:),2);
    VY = flip(vy_corrected(:,x1:x2,:),2);
    v_y= flip(vy_corrected(:,:,:),2);
    VZ = flip(vz_corrected(:,x1:x2,:),2);
    v_z= flip(vz_corrected(:,:,:),2);
    M  = flip(mag_corrected(:,x1:x2,:),2);
    m  = flip(mag_corrected(:,:,:),2);
        
%montage of velocities

lbl_space=zeros(size(VX,1),size(VX,2));

I_montage=cat(3,cat(3,lbl_space,VX),cat(3,lbl_space,VY),cat(3,lbl_space,VZ));
cmap=jet(201);
cmap(101,:) = 0;
border_w=20;
montage(I_montage,'DisplayRange',[-2.5 2.5],'Size',[3,size(VX,3)+1],'BorderSize',[border_w,0],'ThumbnailSize',[size(VX,1),size(VX,2)]);
colormap(cmap)
set(gcf,'color','black');
h=colorbar;
h.Color='w';
h.FontSize=12;
ylabel(h, '[cm/s]','Interpreter','latex','FontSize',20);
set(h,'TickLabelInterpreter','latex','FontSize',20); 

text(round(size(VX,2)/3),round(border_w+size(VX,1)/2),'$v_{x}$','FontSize',32,'Interpreter','latex','Color','white')
text(round(size(VX,2)/3),round(3*border_w+size(VX,1)*1.5),'$v_{y}$','FontSize',32,'Interpreter','latex','Color','white')
text(round(size(VX,2)/3),round(5*border_w+size(VX,1)*2.5),'$v_{z}$','FontSize',32,'Interpreter','latex','Color','white')

filename=['vel_',name_suffix,'.png'];
export_fig(filename,'-m8')

close all

%montage of magnitude

    %intensity correction
    for i=1:size(M,3)
        M(:,:,i)=mat2gray(M(:,:,i));
    end

montage(M,'DisplayRange',[0 1],'Size',[1,size(M,3)],'BorderSize',[border_w,0],'ThumbnailSize',[size(M,1),size(M,2)]);
set(gcf,'color','black');

filename=['mag_',name_suffix,'.png'];
export_fig(filename,'-m8')

close all


end