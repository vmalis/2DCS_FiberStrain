
function struct=strain_tensor_calcs(vx,VX,vy,VY,vz,VZ,xs,ys,zs,pixelspacing,trigger)

clc

%% displacement calculations Euler way

dx=zeros(size(vx));
dy=zeros(size(vy));
dz=zeros(size(vz));

for t=1:size(trigger,1)
    for x=1:size(vx,2)
        for y=1:size(vx,1)
            for z=1:size(vx,3)
                dx(y,x,z,t+1) = dx(y,x,z,t)+VX(y,x,z,t)*trigger(t);
                dy(y,x,z,t+1) = dy(y,x,z,t)+(-1)*VY(y,x,z,t)*trigger(t);
                dz(y,x,z,t+1) = dz(y,x,z,t)+VZ(y,x,z,t)*trigger(t);
            end
        end
     end
end                         
        



%% displacement calculations Lagrange way

dxL=zeros(size(vx));
dyL=zeros(size(vy));
dzL=zeros(size(vz));

for t=1:size(trigger,1)
    for x=1:size(vx,2)
        for y=1:size(vx,1)
            for z=1:size(vx,3)
                
                X=round(xs(y,x,z,t));
                Y=round(ys(y,x,z,t));
                Z=round(zs(y,x,z,t));
                
                % take care of boundary
                if X>size(vx,2)
                    X=size(vx,2);
                end
                
                if Y>size(vx,1)
                    Y=size(vx,1);
                end
                
                if Z>size(vx,3)
                    Z=size(vx,3);
                end
                
                if X<1
                    X=1;
                end
                
                if Y<1
                    Y=1;
                end
                
                if Z<1
                    Z=1;
                end
                
                
                
                dxL(y,x,z,t+1) = dx(y,x,z,t)+vx(Y,X,Z,t)*trigger(t);
                dyL(y,x,z,t+1) = dy(y,x,z,t)+(-1)*vy(Y,X,Z,t)*trigger(t);
                dzL(y,x,z,t+1) = dz(y,x,z,t)+vz(Y,X,Z,t)*trigger(t);
                
            end
        end
     end
end     



% allocate for tensors
E=zeros([size(vx),3,3]);
L=E;
SR=E;
SR_L=E;


%% strain Tensor euler way
%multiplication coefficient to keep same convention as previous papers
k=10;

disp("Strain tensor calculations Euler way")
disp("=======================================")
tic
for i=1:size(trigger,1)+1
    E(:,:,:,i,:,:)=k*tensor_calcs(squeeze(dx(:,:,:,i)),...
                                    squeeze(dy(:,:,:,i)),...
                                            squeeze(dz(:,:,:,i)),...
                                                 pixelspacing);
end

E(1,:,:,:,:,:)=[];
E(end,:,:,:,:,:)=[];
E(:,1,:,:,:,:)=[];
E(:,end,:,:,:,:)=[];
E(:,:,1,:,:,:)=[];
E(:,:,end,:,:,:)=[];
toc


%% strain Tensor Lagrange way
%multiplication coefficient to keep same convention as previous papers
k=10;

disp("Strain tensor calculations Lagrange way")
disp("=======================================")
tic
for i=1:size(trigger,1)+1
    L(:,:,:,i,:,:)=k*tensor_calcs(squeeze(dxL(:,:,:,i)),...
                                    squeeze(dyL(:,:,:,i)),...
                                            squeeze(dzL(:,:,:,i)),...
                                                 pixelspacing);
end

L(1,:,:,:,:,:)=[];
L(end,:,:,:,:,:)=[];
L(:,1,:,:,:,:)=[];
L(:,end,:,:,:,:)=[];
L(:,:,1,:,:,:)=[];
L(:,:,end,:,:,:)=[];
toc




%% strain rate Tensor
k=10000;

disp("Strain rate tensor calculations")
disp("=======================================")
tic
for i=1:size(trigger,1)+1
    SR(:,:,:,i,:,:)=k*tensor_calcs(squeeze(VX(:,:,:,i)),...
                                    squeeze(VY(:,:,:,i)),...
                                            squeeze(VZ(:,:,:,i)),...
                                                 pixelspacing);
end

SR(1,:,:,:,:,:)=[];
SR(end,:,:,:,:,:)=[];
SR(:,1,:,:,:,:)=[];
SR(:,end,:,:,:,:)=[];
SR(:,:,1,:,:,:)=[];
SR(:,:,end,:,:,:)=[];

toc

%% eigenvalue decomposition for 
disp("Strain rate and Strain ev decomposition")
disp("=======================================")
tic


sizeE=size(E);
VE=zeros(sizeE);
VL=zeros(sizeE);
VSR=zeros(sizeE);
LE=zeros(sizeE(1:end-1));
LL=zeros(sizeE(1:end-1));
LSR=zeros(sizeE(1:end-1));
LE_temp=zeros(3);
LL_temp=zeros(3);
LSR_temp=zeros(3);

for x=1:sizeE(1)
    for y=1:sizeE(2)
        for z=1:sizeE(3)
            for frame=1:sizeE(4)
                
             [VE(x,y,z,frame,:,:),LE_temp]=eig(squeeze(E(x,y,z,frame,:,:)));
             LE(x,y,z,frame,:)=[LE_temp(1,1),LE_temp(2,2),LE_temp(3,3)];
             
             [VL(x,y,z,frame,:,:),LL_temp]=eig(squeeze(L(x,y,z,frame,:,:)));
             LL(x,y,z,frame,:)=[LL_temp(1,1),LL_temp(2,2),LL_temp(3,3)];
             
             [VSR(x,y,z,frame,:,:),LSR_temp]=eig(squeeze(SR(x,y,z,frame,:,:)));
             LSR(x,y,z,frame,:)=[LSR_temp(1,1),LSR_temp(2,2),LSR_temp(3,3)];
                
            end
        end
    end
end
toc

%% max shear calaculations (i.e. octahedral shear strain) Euler Way
disp("Shear Strain and Strain rate")
disp("=======================================")
tic

maxShearE=zeros(sizeE(1:4));
maxShearL=zeros(sizeE(1:4));
maxShearSR=zeros(sizeE(1:4));

maxShearE(:,:,:,:)=2/3*sqrt((E(:,:,:,:,1,1)-E(:,:,:,:,2,2)).^2+...
                            +(E(:,:,:,:,1,1)-E(:,:,:,:,3,3)).^2+...
                            +(E(:,:,:,:,2,2)-E(:,:,:,:,3,3)).^2+...
            +6*(E(:,:,:,:,1,2).^2+E(:,:,:,:,1,3).^2+E(:,:,:,:,2,3).^2));
        
        
maxShearL(:,:,:,:)=2/3*sqrt((L(:,:,:,:,1,1)-L(:,:,:,:,2,2)).^2+...
                            +(L(:,:,:,:,1,1)-L(:,:,:,:,3,3)).^2+...
                            +(L(:,:,:,:,2,2)-L(:,:,:,:,3,3)).^2+...
            +6*(L(:,:,:,:,1,2).^2+L(:,:,:,:,1,3).^2+L(:,:,:,:,2,3).^2));        
        
         
maxShearSR(:,:,:,:)=2/3*sqrt((SR(:,:,:,:,1,1)-SR(:,:,:,:,2,2)).^2+...
                            +(SR(:,:,:,:,1,1)-SR(:,:,:,:,3,3)).^2+...
                            +(SR(:,:,:,:,2,2)-SR(:,:,:,:,3,3)).^2+...
            +6*(SR(:,:,:,:,1,2).^2+SR(:,:,:,:,1,3).^2+SR(:,:,:,:,2,3).^2));
toc
%%

%Strain tensor Euler way
struct.E=padarray(E,[1 1],0,'both');

%Strain tensor Lagrange way
struct.L=padarray(L,[1 1],0,'both');

%Strain rate tensor
struct.SR=padarray(SR,[1 1],0,'both');

% d* displacements
struct.dx=dx;
struct.dy=dy;
struct.dz=dz;

% *s tracked coordinates at each frame
struct.xs=xs;
struct.ys=ys;
struct.zs=zs;

% v* velocity along the * direction at each frame
struct.vx=vx;
struct.vy=vy;
struct.vz=vz;

% Strain eigenvalues
struct.E_lambda=padarray(LE,[1 1],0,'both');
% Strain eigen vectors (packed column by column as 3x3 matrix)
struct.E_vector=padarray(VE,[1 1],0,'both');
% max shear strain per voxel in roi
struct.ShearE_max=padarray(maxShearE,[1 1],0,'both');
%volumetric strain i.e. sum of diagonal components   
struct.E_Volumetric=padarray(sum(LE,5),[1 1],0,'both');


% Strain eigenvalues
struct.L_lambda=padarray(LL,[1 1],0,'both');
% Strain eigen vectors (packed column by column as 3x3 matrix)
struct.L_vector=padarray(VL,[1 1],0,'both');
% max shear strain per voxel in roi
struct.ShearL_max=padarray(maxShearL,[1 1],0,'both');
%volumetric strain i.e. sum of diagonal components   
struct.L_Volumetric=padarray(sum(LL,5),[1 1],0,'both');

    
% Strain rate eigenvalues
struct.SR_lambda=padarray(LSR,[1 1],0,'both');
% Strain rate eigen vectors
struct.SR_vector=padarray(VSR,[1 1],0,'both');
% max shear strain rate per voxel in roi
struct.ShearSR_max=padarray(maxShearSR,[1 1],0,'both');
    
end