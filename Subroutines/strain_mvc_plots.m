function []=strain_mvc_plots(struct,roi_name)

dt=0.14;
mkdir(roi_name)
cd(roi_name)


% displacement
limits=[-15, 15];
v=cat(2,struct(1).dx, struct(2).dx, struct(3).dx);
sd=cat(2,struct(1).dx_sd, struct(2).dx_sd, struct(3).dx_sd);
strain_plot(10*v,10*sd,'$\Delta_{x} \quad [\mathrm{mm}]$','dx',limits,dt);

v=cat(2,struct(1).dy, struct(2).dy, struct(3).dy);
sd=cat(2,struct(1).dy_sd, struct(2).dy_sd, struct(3).dy_sd);
strain_plot(10*v,10*sd,'$\Delta_{y} \quad [\mathrm{mm}]$','dy',limits,dt);

v=cat(2,struct(1).dz, struct(2).dz, struct(3).dz);
sd=cat(2,struct(1).dz_sd, struct(2).dz_sd, struct(3).dz_sd);
strain_plot(10*v,10*sd,'$\Delta_{z} \quad [\mathrm{mm}]$','dz',limits,dt);


%velocity
limits=[-5, 5];
v=cat(2,struct(1).vx, struct(2).vx, struct(3).vx);
sd=cat(2,struct(1).vx_sd, struct(2).vx_sd, struct(3).vx_sd);
strain_plot(v,sd,'$v_{x} \quad [\mathrm{cm/s}]$','vx',limits,dt);

v=cat(2,struct(1).vy, struct(2).vy, struct(3).vy);
sd=cat(2,struct(1).vy_sd, struct(2).vy_sd, struct(3).vy_sd);
strain_plot(v,sd,'$v_{y} \quad [\mathrm{cm/s}]$','vy',limits,dt);

v=cat(2,struct(1).vz, struct(2).vz, struct(3).vz);
sd=cat(2,struct(1).vz_sd, struct(2).vz_sd, struct(3).vz_sd);
strain_plot(v,sd,'$v_{z} \quad [\mathrm{cm/s}]$','vz',limits,dt);


% euler strain
limits=[-1, 1];
v=cat(2,struct(1).E_lambda(:,1), struct(2).E_lambda(:,1), struct(3).E_lambda(:,1));
sd=cat(2,struct(1).E_lambda_sd(:,1), struct(2).E_lambda_sd(:,1), struct(3).E_lambda_sd(:,1));
strain_plot(v,sd,'$E_{\lambda_{1}} \quad [\mathrm{mm/mm}]$','E1',limits,dt);

v=cat(2,struct(1).E_lambda(:,2), struct(2).E_lambda(:,2), struct(3).E_lambda(:,2));
sd=cat(2,struct(1).E_lambda_sd(:,2), struct(2).E_lambda_sd(:,2), struct(3).E_lambda_sd(:,2));
strain_plot(v,sd,'$E_{\lambda_{2}} \quad [\mathrm{mm/mm}]$','E2',limits,dt);

v=cat(2,struct(1).E_lambda(:,3), struct(2).E_lambda(:,3), struct(3).E_lambda(:,3));
sd=cat(2,struct(1).E_lambda_sd(:,3), struct(2).E_lambda_sd(:,3), struct(3).E_lambda_sd(:,3));
strain_plot(v,sd,'$E_{\lambda_{3}} \quad [\mathrm{mm/mm}]$','E3',limits,dt);

v=cat(2,struct(1).ShearE_max, struct(2).ShearE_max, struct(3).ShearE_max);
sd=cat(2,struct(1).ShearE_max_sd(:,1), struct(2).ShearE_max_sd(:,1), struct(3).ShearE_max_sd(:,1));
strain_plot(v,sd,'$E_{max} \quad [\mathrm{mm/mm}]$','Emax',limits,dt);

v=cat(2,struct(1).E_Volumetric, struct(2).E_Volumetric, struct(3).E_Volumetric);
sd=cat(2,struct(1).E_Volumetric_sd, struct(2).E_Volumetric_sd, struct(3).E_Volumetric_sd);
strain_plot(v,sd,'$E_{vol} \quad [\mathrm{mm^3/mm^3}]$','EV',limits,dt);

% lagrangian strain
limits=[-1, 1];
v=cat(2,struct(1).L_lambda(:,1), struct(2).L_lambda(:,1), struct(3).L_lambda(:,1));
sd=cat(2,struct(1).L_lambda_sd(:,1), struct(2).L_lambda_sd(:,1), struct(3).L_lambda_sd(:,1));
strain_plot(v,sd,'$L_{\lambda_{1}} \quad [\mathrm{mm/mm}]$','L1',limits,dt);

v=cat(2,struct(1).L_lambda(:,2), struct(2).L_lambda(:,2), struct(3).L_lambda(:,2));
sd=cat(2,struct(1).L_lambda_sd(:,2), struct(2).L_lambda_sd(:,2), struct(3).L_lambda_sd(:,2));
strain_plot(v,sd,'$L_{\lambda_{2}} \quad [\mathrm{mm/mm}]$','L2',limits,dt);

v=cat(2,struct(1).L_lambda(:,3), struct(2).L_lambda(:,3), struct(3).L_lambda(:,3));
sd=cat(2,struct(1).L_lambda_sd(:,3), struct(2).L_lambda_sd(:,3), struct(3).L_lambda_sd(:,3));
strain_plot(v,sd,'$L_{\lambda_{3}} \quad [\mathrm{mm/mm}]$','L3',limits,dt);

v=cat(2,struct(1).ShearL_max, struct(2).ShearL_max, struct(3).ShearL_max);
sd=cat(2,struct(1).ShearL_max_sd(:,1), struct(2).ShearL_max_sd(:,1), struct(3).ShearL_max_sd(:,1));
strain_plot(v,sd,'$L_{max} \quad [\mathrm{mm/mm}]$','Lmax',limits,dt);

v=cat(2,struct(1).L_Volumetric, struct(2).L_Volumetric, struct(3).L_Volumetric);
sd=cat(2,struct(1).L_Volumetric_sd, struct(2).L_Volumetric_sd, struct(3).L_Volumetric_sd);
strain_plot(v,sd,'$L_{vol} \quad [\mathrm{mm^3/mm^3}]$','LV',limits,dt);

% euler strain rate
limits=[-2000,2000];
v=cat(2,struct(1).SR_E_lambda(:,1), struct(2).SR_E_lambda(:,1), struct(3).SR_E_lambda(:,1));
sd=cat(2,struct(1).SR_E_lambda_sd(:,1), struct(2).SR_E_lambda_sd(:,1), struct(3).SR_E_lambda_sd(:,1));
strain_plot(v,sd,'$SR_{\lambda_{1}} \quad [\mathrm{m s^{-1}}]$','SRE1',limits,dt);

v=cat(2,struct(1).SR_E_lambda(:,2), struct(2).SR_E_lambda(:,2), struct(3).SR_E_lambda(:,2));
sd=cat(2,struct(1).SR_E_lambda_sd(:,2), struct(2).SR_E_lambda_sd(:,2), struct(3).SR_E_lambda_sd(:,2));
strain_plot(v,sd,'$SR_{\lambda_{2}} \quad [\mathrm{m s^{-1}}]$','SRE2',limits,dt);

v=cat(2,struct(1).SR_E_lambda(:,3), struct(2).SR_E_lambda(:,3), struct(3).SR_E_lambda(:,3));
sd=cat(2,struct(1).SR_E_lambda_sd(:,3), struct(2).SR_E_lambda_sd(:,3), struct(3).SR_E_lambda_sd(:,3));
strain_plot(v,sd,'$SR_{\lambda_{3}} \quad [\mathrm{m s^{-1}}]$','SRE3',limits,dt);

limits=[-2000, 2000];
v=cat(2,struct(1).ShearSR_E_max, struct(2).ShearSR_E_max, struct(3).ShearSR_E_max);
sd=cat(2,struct(1).ShearSR_E_max_sd(:,1), struct(2).ShearSR_E_max_sd(:,1), struct(3).ShearSR_E_max_sd(:,1));
strain_plot(v,sd,'$SR_{max} \quad [\mathrm{m s^{-1}}]$','SREmax',limits,dt);


% lagrangian strain rate Lagrangian
limits=[-2000,2000];

v=cat(2,struct(1).SR_lambda(:,1), struct(2).SR_lambda(:,1), struct(3).SR_lambda(:,1));
sd=cat(2,struct(1).SR_lambda_sd(:,1), struct(2).SR_lambda_sd(:,1), struct(3).SR_lambda_sd(:,1));
strain_plot(v,sd,'$SR_{\lambda_{1}} \quad [\mathrm{m s^{-1}}]$','SR1',limits,dt);

v=cat(2,struct(1).SR_lambda(:,2), struct(2).SR_lambda(:,2), struct(3).SR_lambda(:,2));
sd=cat(2,struct(1).SR_lambda_sd(:,2), struct(2).SR_lambda_sd(:,2), struct(3).SR_lambda_sd(:,2));
strain_plot(v,sd,'$SR_{\lambda_{2}} \quad [\mathrm{m s^{-1}}]$','SR2',limits,dt);

v=cat(2,struct(1).SR_lambda(:,3), struct(2).SR_lambda(:,3), struct(3).SR_lambda(:,3));
sd=cat(2,struct(1).SR_lambda_sd(:,3), struct(2).SR_lambda_sd(:,3), struct(3).SR_lambda_sd(:,3));
strain_plot(v,sd,'$SR_{\lambda_{3}} \quad [\mathrm{m s^{-1}}]$','SR3',limits,dt);


limits=[-3000, 3000];
v=cat(2,struct(1).ShearSR_max, struct(2).ShearSR_max, struct(3).ShearSR_max);
sd=cat(2,struct(1).ShearSR_max_sd(:,1), struct(2).ShearSR_max_sd(:,1), struct(3).ShearSR_max_sd(:,1));
strain_plot(v,sd,'$SR_{max} \quad [\mathrm{m s^{-1}}]$','SRmax',limits,dt);



%%% fiber aligned

% euler strain
limits=[-1, 1];

v=cat(2,struct(1).ED(:,1,1), struct(2).ED(:,1,1), struct(3).ED(:,1,1));
sd=cat(2,struct(1).ED_sd(:,1,1), struct(2).ED_sd(:,1,1), struct(3).ED_sd(:,1,1));
strain_plot(v,sd,'$E_{ff} \quad [\mathrm{mm/mm}]$','Eff',limits,dt);

v=cat(2,struct(1).ED(:,2,2), struct(2).ED(:,2,2), struct(3).ED(:,2,2));
sd=cat(2,struct(1).ED_sd(:,1,1), struct(2).ED_sd(:,1,1), struct(3).ED_sd(:,1,1));
strain_plot(v,sd,'$E_{ss} \quad [\mathrm{mm/mm}]$','Ess',limits,dt);

v=cat(2,struct(1).ED(:,3,3), struct(2).ED(:,3,3), struct(3).ED(:,3,3));
sd=cat(2,struct(1).ED_sd(:,3,3), struct(2).ED_sd(:,3,3), struct(3).ED_sd(:,3,3));
strain_plot(v,sd,'$E_{tt} \quad [\mathrm{mm/mm}]$','Ett',limits,dt);

v=cat(2,struct(1).ED(:,1,2), struct(2).ED(:,1,2), struct(3).ED(:,1,2));
sd=cat(2,struct(1).ED_sd(:,1,2), struct(2).ED_sd(:,1,2), struct(3).ED_sd(:,1,2));
strain_plot(v,sd,'$E_{fs} \quad [\mathrm{mm/mm}]$','Efs',limits,dt);

v=cat(2,struct(1).ED(:,1,3), struct(2).ED(:,1,3), struct(3).ED(:,1,3));
sd=cat(2,struct(1).ED_sd(:,1,3), struct(2).ED_sd(:,1,3), struct(3).ED_sd(:,1,3));
strain_plot(v,sd,'$E_{ft} \quad [\mathrm{mm/mm}]$','Eft',limits,dt);

v=cat(2,struct(1).ED(:,2,3), struct(2).ED(:,2,3), struct(3).ED(:,2,3));
sd=cat(2,struct(1).ED_sd(:,2,3), struct(2).ED_sd(:,2,3), struct(3).ED_sd(:,2,3));
strain_plot(v,sd,'$E_{st} \quad [\mathrm{mm/mm}]$','Est',limits,dt);

% lagrangian strain
limits=[-1, 1];

v=cat(2,struct(1).LD(:,1,1), struct(2).LD(:,1,1), struct(3).LD(:,1,1));
sd=cat(2,struct(1).LD_sd(:,1,1), struct(2).LD_sd(:,1,1), struct(3).LD_sd(:,1,1));
strain_plot(v,sd,'$L_{ff} \quad [\mathrm{mm/mm}]$','Lff',limits,dt);

v=cat(2,struct(1).LD(:,2,2), struct(2).LD(:,2,2), struct(3).LD(:,2,2));
sd=cat(2,struct(1).LD_sd(:,1,1), struct(2).LD_sd(:,1,1), struct(3).LD_sd(:,1,1));
strain_plot(v,sd,'$L_{ss} \quad [\mathrm{mm/mm}]$','Lss',limits,dt);

v=cat(2,struct(1).LD(:,3,3), struct(2).LD(:,3,3), struct(3).LD(:,3,3));
sd=cat(2,struct(1).LD_sd(:,3,3), struct(2).LD_sd(:,3,3), struct(3).LD_sd(:,3,3));
strain_plot(v,sd,'$L_{tt} \quad [\mathrm{mm/mm}]$','Ltt',limits,dt);

v=cat(2,struct(1).LD(:,1,2), struct(2).LD(:,1,2), struct(3).LD(:,1,2));
sd=cat(2,struct(1).LD_sd(:,1,2), struct(2).LD_sd(:,1,2), struct(3).LD_sd(:,1,2));
strain_plot(v,sd,'$L_{fs} \quad [\mathrm{mm/mm}]$','Lfs',limits,dt);

v=cat(2,struct(1).LD(:,1,3), struct(2).LD(:,1,3), struct(3).LD(:,1,3));
sd=cat(2,struct(1).LD_sd(:,1,3), struct(2).LD_sd(:,1,3), struct(3).LD_sd(:,1,3));
strain_plot(v,sd,'$L_{ft} \quad [\mathrm{mm/mm}]$','Lft',limits,dt);

v=cat(2,struct(1).LD(:,2,3), struct(2).LD(:,2,3), struct(3).LD(:,2,3));
sd=cat(2,struct(1).LD_sd(:,2,3), struct(2).LD_sd(:,2,3), struct(3).LD_sd(:,2,3));
strain_plot(v,sd,'$L_{st} \quad [\mathrm{mm/mm}]$','Lst',limits,dt);


% lagrangian strain rate
limits=[-2000, 2000];
v=cat(2,struct(1).SRD(:,1,1), struct(2).SRD(:,1,1), struct(3).SRD(:,1,1));
sd=cat(2,struct(1).SRD_sd(:,1,1), struct(2).SRD_sd(:,1,1), struct(3).SRD_sd(:,1,1));
strain_plot(v,sd,'$SR_{ff} \quad [\mathrm{m s^{-1}}]$','SRff',limits,dt);

v=cat(2,struct(1).SRD(:,2,2), struct(2).SRD(:,2,2), struct(3).SRD(:,2,2));
sd=cat(2,struct(1).SRD_sd(:,1,1), struct(2).SRD_sd(:,1,1), struct(3).SRD_sd(:,1,1));
strain_plot(v,sd,'$SR_{ss} \quad [\mathrm{m s^{-1}}]$','SRss',limits,dt);

v=cat(2,struct(1).SRD(:,3,3), struct(2).SRD(:,3,3), struct(3).SRD(:,3,3));
sd=cat(2,struct(1).SRD_sd(:,3,3), struct(2).SRD_sd(:,3,3), struct(3).SRD_sd(:,3,3));
strain_plot(v,sd,'$SR_{tt} \quad [\mathrm{m s^{-1}}]$','SRtt',limits,dt);

v=cat(2,struct(1).SRD(:,1,2), struct(2).SRD(:,1,2), struct(3).SRD(:,1,2));
sd=cat(2,struct(1).SRD_sd(:,1,2), struct(2).SRD_sd(:,1,2), struct(3).SRD_sd(:,1,2));
strain_plot(v,sd,'$SR_{fs} \quad [\mathrm{m s^{-1}}]$','SRfs',limits,dt);

v=cat(2,struct(1).SRD(:,1,3), struct(2).SRD(:,1,3), struct(3).SRD(:,1,3));
sd=cat(2,struct(1).SRD_sd(:,1,3), struct(2).SRD_sd(:,1,3), struct(3).SRD_sd(:,1,3));
strain_plot(v,sd,'$SR_{ft} \quad [\mathrm{m s^{-1}}]$','SRft',limits,dt);

v=cat(2,struct(1).SRD(:,2,3), struct(2).SRD(:,2,3), struct(3).SRD(:,2,3));
sd=cat(2,struct(1).SRD_sd(:,2,3), struct(2).SRD_sd(:,2,3), struct(3).SRD_sd(:,2,3));
strain_plot(v,sd,'$SR_{st} \quad [\mathrm{m s^{-1}}]$','SRst',limits,dt);


% euler strain rate
limits=[-2000, 2000];
v=cat(2,struct(1).SRD_E(:,1,1), struct(2).SRD_E(:,1,1), struct(3).SRD_E(:,1,1));
sd=cat(2,struct(1).SRD_E_sd(:,1,1), struct(2).SRD_E_sd(:,1,1), struct(3).SRD_E_sd(:,1,1));
strain_plot(v,sd,'$SR_{ff} \quad [\mathrm{m s^{-1}}]$','SRff',limits,dt);

v=cat(2,struct(1).SRD_E(:,2,2), struct(2).SRD_E(:,2,2), struct(3).SRD_E(:,2,2));
sd=cat(2,struct(1).SRD_E_sd(:,1,1), struct(2).SRD_E_sd(:,1,1), struct(3).SRD_E_sd(:,1,1));
strain_plot(v,sd,'$SR_{ss} \quad [\mathrm{m s^{-1}}]$','SRss',limits,dt);

v=cat(2,struct(1).SRD_E(:,3,3), struct(2).SRD_E(:,3,3), struct(3).SRD_E(:,3,3));
sd=cat(2,struct(1).SRD_E_sd(:,3,3), struct(2).SRD_E_sd(:,3,3), struct(3).SRD_E_sd(:,3,3));
strain_plot(v,sd,'$SR_{tt} \quad [\mathrm{m s^{-1}}]$','SRtt',limits,dt);

v=cat(2,struct(1).SRD_E(:,1,2), struct(2).SRD_E(:,1,2), struct(3).SRD_E(:,1,2));
sd=cat(2,struct(1).SRD_E_sd(:,1,2), struct(2).SRD_E_sd(:,1,2), struct(3).SRD_E_sd(:,1,2));
strain_plot(v,sd,'$SR_{fs} \quad [\mathrm{m s^{-1}}]$','SRfs',limits,dt);

v=cat(2,struct(1).SRD_E(:,1,3), struct(2).SRD_E(:,1,3), struct(3).SRD_E(:,1,3));
sd=cat(2,struct(1).SRD_E_sd(:,1,3), struct(2).SRD_E_sd(:,1,3), struct(3).SRD_E_sd(:,1,3));
strain_plot(v,sd,'$SR_{ft} \quad [\mathrm{m s^{-1}}]$','SRft',limits,dt);

v=cat(2,struct(1).SRD_E(:,2,3), struct(2).SRD_E(:,2,3), struct(3).SRD_E(:,2,3));
sd=cat(2,struct(1).SRD_E_sd(:,2,3), struct(2).SRD_E_sd(:,2,3), struct(3).SRD_E_sd(:,2,3));
strain_plot(v,sd,'$SR_{st} \quad [\mathrm{m s^{-1}}]$','SRst',limits,dt);

cd ..

end

