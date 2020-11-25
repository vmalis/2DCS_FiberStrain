function struct = strain_ROI(structure,mask)


mask(mask==0)=NaN;
n_frames=size(mask,4);

mask0                   = repmat(mask(:,:,1),[1,1,n_frames]);
mask_vector_0           = repmat(mask0,[1,1,1,3]);
mask_vector_dynamic     = repmat(mask,[1,1,1,1,3]);
mask_tensor_0           = repmat(mask0,[1,1,1,3,3]);
mask_tensor_dynamic     = repmat(mask,[1,1,1,1,3,3]);

%displacements
temp_dx                 =  squeeze(structure.dx(:,:,2,:)).*mask0;
struct.dx               =  squeeze(nanmean(nanmean(temp_dx,1),2));
struct.dx_sd            =  squeeze(nanstd(nanstd(temp_dx,1,1),1,2));

temp_dy                 =  squeeze(structure.dy(:,:,2,:)).*mask0;
struct.dy               =  squeeze(nanmean(nanmean(temp_dy,1),2));
struct.dy_sd            =  squeeze(nanstd(nanstd(temp_dy,1,1),1,2));

temp_dz                 =  squeeze(structure.dz(:,:,2,:)).*mask0;
struct.dz               =  squeeze(nanmean(nanmean(temp_dz,1),2));
struct.dz_sd            =  squeeze(nanstd(nanstd(temp_dz,1,1),1,2));

temp_dr                 =  squeeze(sqrt(temp_dx.^2 + temp_dy.^2 + temp_dz.^2));
struct.dr               =  squeeze(nanmean(nanmean(temp_dr,1),2));
struct.dr_sd            =  squeeze(nanstd(nanstd(temp_dr,1,1),1,2));


%velocities
temp_vx                 =  squeeze(structure.vx(:,:,2,:)).*mask0;
struct.vx               =  squeeze(nanmean(nanmean(temp_vx,1),2));
struct.vx_sd            =  squeeze(nanstd(nanstd(temp_vx,1,1),1,2));

temp_vy                 =  squeeze(structure.vy(:,:,2,:)).*mask0;
struct.vy               =  squeeze(nanmean(nanmean(temp_vy,1),2));
struct.vy_sd            =  squeeze(nanstd(nanstd(temp_vy,1,1),1,2));

temp_vz                 =  squeeze(structure.vz(:,:,2,:)).*mask0;
struct.vz               =  squeeze(nanmean(nanmean(temp_vz,1),2));
struct.vz_sd            =  squeeze(nanstd(nanstd(temp_vz,1,1),1,2));

temp_vr                 =  squeeze(sqrt(temp_vx.^2 + temp_vy.^2 + temp_vz.^2));
struct.vr               =  squeeze(nanmean(nanmean(temp_vr,1),2));
struct.vr_sd            =  squeeze(nanstd(nanstd(temp_vr,1,1),1,2));



%strain Euler way
temp_E_lambda           = squeeze(structure.E_lambda).*mask_vector_0;
struct.E_lambda         = squeeze(nanmean(nanmean(temp_E_lambda,2),1));
struct.E_lambda_sd      = squeeze(nanstd(nanstd(temp_E_lambda,1,1),1,2));

temp_ShearE_max         = squeeze(structure.ShearE_max).*mask0;
struct.ShearE_max       = squeeze(nanmean(nanmean(temp_ShearE_max,2),1));
struct.ShearE_max_sd    = squeeze(nanstd(nanstd(temp_ShearE_max,1,1),1,2));

temp_E_Volumetric       = squeeze(structure.E_Volumetric).*mask0;
struct.E_Volumetric     = squeeze(nanmean(nanmean(temp_E_Volumetric,2),1));
struct.E_Volumetric_sd  = squeeze(nanstd(nanstd(temp_E_Volumetric,1,1),1,2));


%strain Lagrange way
temp_L_lambda           = squeeze(structure.L_lambda).*mask_vector_0;
struct.L_lambda         = squeeze(nanmean(nanmean(temp_L_lambda,2),1));
struct.L_lambda_sd      = squeeze(nanstd(nanstd(temp_L_lambda,1,1),1,2));

temp_ShearL_max         = squeeze(structure.ShearL_max).*mask0;
struct.ShearL_max       = squeeze(nanmean(nanmean(temp_ShearL_max,2),1));
struct.ShearL_max_sd    = squeeze(nanstd(nanstd(temp_ShearL_max,1,1),1,2));

temp_L_Volumetric       = squeeze(structure.L_Volumetric).*mask0;
struct.L_Volumetric     = squeeze(nanmean(nanmean(temp_L_Volumetric,2),1));
struct.L_Volumetric_sd  = squeeze(nanstd(nanstd(temp_L_Volumetric,1,1),1,2));

%strain rate Lagrange way, our regular way

temp_SR_lambda          = squeeze(structure.SR_lambda).*squeeze(mask_vector_dynamic(:,:,2,:,:));
struct.SR_lambda        = squeeze(nanmean(nanmean(temp_SR_lambda,2),1));
struct.SR_lambda_sd     = squeeze(nanstd(nanstd(temp_SR_lambda,1,1),1,2));

temp_ShearSR_max        = squeeze(structure.ShearSR_max).*squeeze(mask(:,:,2,:));
struct.ShearSR_max      = squeeze(nanmean(nanmean(temp_ShearSR_max,2),1));
struct.ShearSR_max_sd   = squeeze(nanstd(nanstd(temp_ShearSR_max,1,1),1,2));


%strain rate Euler way (static mask)
temp_SR_E_lambda        = squeeze(structure.SR_lambda).*mask_vector_0;
struct.SR_E_lambda      = squeeze(nanmean(nanmean(temp_SR_E_lambda,2),1));
struct.SR_E_lambda_sd   = squeeze(nanstd(nanstd(temp_SR_E_lambda,1,1),1,2));

temp_ShearSR_E_max      = squeeze(structure.ShearSR_max).*mask0;
struct.ShearSR_E_max    = squeeze(nanmean(nanmean(temp_ShearSR_E_max,2),1));
struct.ShearSR_E_max_sd = squeeze(nanstd(nanstd(temp_ShearSR_E_max,1,1),1,2));

%fiber alligned

% strain Euler way
temp_ED                 = squeeze(structure.ED).*mask_tensor_0;
struct.ED               = squeeze(nanmean(nanmean(temp_ED,2),1));
struct.ED_sd            = squeeze(nanstd(nanstd(temp_ED,1,1),1,2));


% strain Lagrange way
temp_LD                 = squeeze(structure.LD).*mask_tensor_0;
struct.LD               = squeeze(nanmean(nanmean(temp_LD,2),1));
struct.LD_sd            = squeeze(nanstd(nanstd(temp_LD,1,1),1,2));

% strain rate Lagrange way

temp_SRD                = squeeze(structure.SRD).*squeeze(mask_tensor_dynamic(:,:,2,:,:,:));
struct.SRD              = squeeze(nanmean(nanmean(temp_SRD,2),1));
struct.SRD_sd           = squeeze(nanstd(nanstd(temp_SRD,1,1),1,2));

% strain rate Euler way
temp_SRD_E              = squeeze(structure.SRD).*mask_tensor_0;
struct.SRD_E            = squeeze(nanmean(nanmean(temp_SRD_E,2),1));
struct.SRD_E_sd         = squeeze(nanstd(nanstd(temp_SRD_E,1,1),1,2));


end