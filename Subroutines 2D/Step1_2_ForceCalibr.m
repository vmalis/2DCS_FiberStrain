%=========================================================================
%   Step1_2 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strane Rate Toolkit
%=========================================================================
%
% This script
%     applies calibration to the force measurments
%
% INput: force_data.mat   
%
% OUTput: force_data.mat (updated with force in Newtons)
%         force plot in Newtons
%__________________________________________________________________________
% required subroutines:
%
% written by Vadim Malis
% 06/15 at UCSD RIL
%==========================================================================
% window_info=sprintf('Choose Data_VEPC folder');
% PathName = uigetdir('~/Desktop',window_info);
% cd(PathName)

file_list=rdir('**/force_data.mat');
list=cell(size(file_list));

for i=1:size(list,1)
  name=file_list(i).name;
  indexes=strfind(file_list(i).name,'patient');
  list{i,1}=name(1:7);%name(5:indexes(1)-2);
end

[uq_list,ia,ic]=unique(list,'stable');
x = {'Purple_main' 'Green_old_big' 'New Green' 'Pedal_old' 'Pedal_new'};
ia=cat(1,ia,size(file_list,1)+1);    

% for i=1:size(uq_list,1)
%     uq_list{i}=[uq_list{i}(end-1:end) ' | ' uq_list{i}(1:end-3)];
% end

% Column names and column format
columnname = {' ','type'};
columnformat = {'char',x};

% Define the data used in uitable; default values are true
d=repmat(x,size(uq_list));
d(:,2:5)=[];
d=cat(2,uq_list,d);
f = figure('CloseRequestFcn',@ediTable);

% Create the uitable
t = uitable('Data', d,... 
            'ColumnName', columnname,...
            'RowName',[],...
            'ColumnFormat', columnformat,...
            'ColumnEditable', [false true],...
            'ColumnWidth',{150,150});

% Set width height and position
ScreenSize=get(0,'Screensize');
tableextent = get(t,'Extent');
oldposition = get(t,'Position');
tnewposition = [oldposition(1) oldposition(2) tableextent(3) tableextent(4)];
set(t, 'Position', tnewposition);
set(f, 'Position',[ScreenSize(3)/2-(tableextent(3)+80)/2 ...
    ScreenSize(4)/2-(tableextent(4)+80)/2 ...
    tableextent(3)+40 tableextent(4)+40]);


set(f,'Tag', 'tableFigure')
set(t,'Tag','table')

% get for table to close and deliver selected data
waitfor(f)


%
bpm=20;

%sample rate
sr=200;

%cycle length in seconds
T=60/bpm;

%create time array for plots
t=0:1/sr:T;
t(end)=[];


selected_calibr=selected_data(:,2);
calibr=zeros(size(file_list));
for i=1:size(uq_list,1)
    
    if strcmp(selected_calibr{i},'Purple_main')        
    calibr(ia(i):ia(i+1)-1)=1;

    elseif strcmp(selected_calibr{i},'Green_old_big')        
    calibr(ia(i):ia(i+1)-1)=2;
    
     elseif strcmp(selected_calibr{i},'New Green')        
    calibr(ia(i):ia(i+1)-1)=3;
    
    elseif strcmp(selected_calibr{i},'Pedal_old')        
    calibr(ia(i):ia(i+1)-1)=4;
    
    elseif strcmp(selected_calibr{i},'Pedal_new')        
    calibr(ia(i):ia(i+1)-1)=5;

    end
end

clearvars -except calibr file_list t T uq_list ia

%-------------applying calibration-----------------------------------------

multiWaitbar('Progress...', 0, 'Color', 'r');
n=size(calibr,1);


force=zeros(600,size(file_list,1));

for i=1:n
    
load(file_list(i).name)

%     if calibr(i)==1
%     k=0.0213;
%     b=0.1329;
%         
%     elseif calibr(i)==2
%     k=0.0089;
%     b=0.065;
%     
%     elseif calibr(i)==3
%     k=0.278;
%     b=-0.0803;
%     
%     
%     elseif calibr(i)==4
%     k=0.0273;
%     b=0.2535;
%     
%     elseif calibr(i)==5
%     k=0.0874;
%     b=-0.0714;


    if calibr(i)==1
    k=46.52;
    b=6.8342;
        
    elseif calibr(i)==2
    k=39.298;
    b=2.811;
    
    elseif calibr(i)==3
    k=111.53;
    b=8.4573;
    
    
    elseif calibr(i)==4
    k=36.5;
    b=9.6101;
    
    elseif calibr(i)==5
    k=35.208;
    b=2.9439;


    end

%     force_mean=(force_mean+b)/k;
%     force_std=(force_std)/k;
%     MAX_force=(MAX_force+b)/k;

force_mean=force_mean*k+b;
force_std=force_std*k+b;
MAX_force=MAX_force*k+b;



force_high=force_mean+force_std;
force_low=force_mean-force_std;    
mvc_line=ones(size(force_mean))*MAX_force;  


%save force plot
h=figure;
set(gcf,'Visible','off');
plot(t,mvc_line,'--m','LineWidth', 1);
hold on
plot(t,force_mean,'-b','LineWidth', 2);
plot(t,force_high,'-r','LineWidth', 1);
plot(t,force_low,'-g','LineWidth', 1);

legend('MVC', 'Mean', 'Upper Bound', 'Lower Bound', 'Location','east')
ylabel('Newtons [N]');
xlabel('Time [s]');
axis([0 T 0-max(force_std) MAX_force*1.1])
axis square
grid on
filename=[file_list(i).name(1:end-14) 'force_average_calibrated.eps'];

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 5 5]);
set(gcf,'color','w');
export_fig(filename)
%print(h,'-dpng',filename);    
    
filename=[file_list(i).name(1:end-14) 'force_data_calibrated.mat']; 
save (filename,'MAX_force')
save (filename,'force_mean','-append')
save (filename,'force_std','-append')
save (filename,'frame','-append')
save (filename,'force_ref','-append')
  
force(:,i)=force_mean;

multiWaitbar('Progress...', i/n);

end

multiWaitbar('Progress...', 'Close');

clearvars -except force ia uq_list