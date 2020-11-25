% Part of ULLS analysis
% minor script to determine the largest crossection of the TS-re muscle

% read all the avaliable processed sets
window_info=sprintf('Choose dti processed data folder');
PathName = uigetdir('~/Desktop',window_info);
cd(PathName)

file_list=rdir('**/Diffusion_data.mat');


%structure to keep crossection information
DATA.ID      = 0;
DATA.time    = 0;
DATA.max_TS  = 0;
DATA.max_GM  = 0;
DATA.max_GL  = 0;
DATA.max_SOL = 0;
DATA.max_TOT = 0;

%pixels to mm^2:  FOV=20cm, pixels resolution = 256
conversion_factor=200*200/256/256;

%processing loop
for i=1:size(file_list)

%fill up ID and time info based on file data    
idx_str=strfind(file_list(i).name,'/');
subject_name=file_list(i).name(1:idx_str(1)-1);
subject_time=file_list(i).name(idx_str(1)+1:idx_str(2)-1);
    
    
DATA(i).ID  = subject_name;
DATA(i).time= subject_time;


load(file_list(i).name,'GM')
load(file_list(i).name,'GL')
load(file_list(i).name,'SOL')

GL(isnan(GL))=0;
GM(isnan(GM))=0;
SOL(isnan(SOL))=0;

MASK=GL+GM+SOL;


[GL_volume,GL_slice]=max(squeeze(squeeze(sum(sum(GL,1),2))));
[GM_volume,GM_slice]=max(squeeze(squeeze(sum(sum(GM,1),2))));
[SOL_volume,SOL_slice]=max(squeeze(squeeze(sum(sum(SOL,1),2))));
[TS_volume,TS_slice]=max(squeeze(squeeze(sum(sum(MASK,1),2))));



DATA(i).max_TS  = TS_volume*conversion_factor;
DATA(i).max_GM  = GM_volume*conversion_factor;
DATA(i).max_GL  = GL_volume*conversion_factor;
DATA(i).max_SOL = SOL_volume*conversion_factor;
DATA(i).max_TOT = DATA(i).max_GM+DATA(i).max_GL+DATA(i).max_SOL;
   
    
end