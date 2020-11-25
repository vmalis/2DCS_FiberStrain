%% -----------------------------------------------------------------------
%  UC San Diego / November 2020 / Vadim Malis
%
%   ISMRM 2021 Strain 2d analysis
%
%-------------------------------------------------------------------------

clc
clear all

disp("Starting reconstruction")
disp("=======================================")

%% read pfile
disp("Step 1: reading pfile to sturcture")
disp("========================================")
raw=pfile2struct;


%% recon

Recon = reconCS(raw(1),'full',1);