function struct = peak_ROI(structure)

% get peak frame
[~,peak_frame]=min(structure.SR_lambda(1:9,1));

% displacement
[struct.dx, struct.dx_sd] =   peak_strains(structure.dx,structure.dx_sd);
[struct.dy, struct.dy_sd] =   peak_strains(structure.dy,structure.dy_sd);
[struct.dz, struct.dz_sd] =   peak_strains(structure.dz,structure.dz_sd);
[struct.dr, struct.dr_sd] =   peak_strains(structure.dr,structure.dr_sd);

% velocities
[struct.vx, struct.vx_sd] =   peak_strains(structure.vx,structure.vx_sd);
[struct.vy, struct.vy_sd] =   peak_strains(structure.vy,structure.vy_sd);
[struct.vz, struct.vz_sd] =   peak_strains(structure.vz,structure.vz_sd);
[struct.vr, struct.vr_sd] =   peak_strains(structure.vr,structure.vr_sd);

struct.vx=abs(struct.vx);
struct.vy=abs(struct.vy);
struct.vz=abs(struct.vz);


%strains euler
[struct.E_lambda, struct.E_lambda_sd] = peak_strains(structure.E_lambda,structure.E_lambda_sd);
[struct.ShearE_max, struct.ShearE_max_sd] = peak_strains(structure.ShearE_max,structure.ShearE_max_sd);
[struct.E_volumetric, struct.E_volumetric_sd] = peak_strains(structure.E_Volumetric,structure.E_Volumetric_sd);

%strains lagrange
[struct.L_lambda, struct.L_lambda_sd] = peak_strains(structure.L_lambda,structure.L_lambda_sd);
[struct.ShearL_max, struct.ShearL_max_sd] = peak_strains(structure.ShearL_max,structure.ShearL_max_sd);
[struct.L_volumetric, struct.L_volumetric_sd] = peak_strains(structure.L_Volumetric,structure.L_Volumetric_sd);

%strain rate lagrange
[struct.SR_lambda, struct.SR_lambda_sd] = peak_strains(structure.SR_lambda,structure.SR_lambda_sd,peak_frame);
[struct.ShearSR_max, struct.ShearSR_max_sd] = peak_strains(structure.ShearSR_max,structure.ShearSR_max_sd,peak_frame);

%strain rate euler
[struct.SR_E_lambda, struct.SR_E_lambda_sd] = peak_strains(structure.SR_lambda,structure.SR_lambda_sd,peak_frame);
[struct.ShearSR_E_max, struct.ShearSR_E_max_sd] = peak_strains(structure.ShearSR_E_max,structure.ShearSR_E_max_sd,peak_frame);


%% fiber aligned

% strain euler
[struct.ED, struct.ED_sd] = peak_strains(structure.ED,structure.ED_sd);

% strain lagrange
[struct.LD, struct.LD_sd] = peak_strains(structure.LD,structure.LD_sd);

% strain rate lagrange
[struct.SRD, struct.SRD_sd] = peak_strains(structure.SRD,structure.SRD_sd,peak_frame);

% strain rate euler
[struct.SRD_E, struct.SRD_E_sd] = peak_strains(structure.SRD_E,structure.SRD_E_sd,peak_frame);



end

