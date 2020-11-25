function [filename] = force_analysis(path,numphases,bpm,permute_flag)

%==========================================================================
%  Subroutine to perform force recording analysis
%==========================================================================
%
% Read *.csv output from Force_Guide.vi 
%
%   data are in units of VOLTS!
%   to get Newtons calibration has to be done
%
% 
% INput:        MVC.csv
%               force_record.csv
%
% OUTput:       *.mat file
%                       1) avaraged force
%                       2) mvc force
%                       3) frame with the peak force
%                       4) avaraged force plot.png saved to the folder
%
%--------------------------------------------------------------------------
% written by Vadim Malis
% 12/14 at UCSD RIL
%==========================================================================

cd(path)

cd ..  % going up to read MVC

%sample rate
sr=200;

%cycle length in seconds
T=60/bpm;

%read folder file structure
MVC_csv=dir('*.csv');

%% MVC
MVC=csvread(MVC_csv.name);
MVC(:,1)=[];
MVC(:,2:end)=[];

%min peak distance (assume we have 3 peaks)
min_peak_distance=round(size(MVC,1)/3);

[MVC_peak]=findpeaks(MVC,'NPeaks',3, 'MinPeakDistance',min_peak_distance);


idx=find(MVC_peak<0.5*max(MVC_peak));
MVC_peak(idx)=[];

MAX_force=mean(MVC_peak);




%% going back to the folder to read force recording
cd(path)
csv_list=dir('*.csv');

%% force data
force_temp=csvread(csv_list(1).name);
trigger=force_temp(:,3);
force=force_temp(:,2);
time=force_temp(:,1);
force=smooth(force,100,'sgolay');
clear force_temp
force=force-min(force);

%save force plot
h=figure;
set(gcf,'Visible','off');
plot(time,force,'-r',time,trigger,'-g');
legend('Force','Trigger');
ylabel('Voltage [V]');
xlabel('Time [s]');

axis([0 max(time) min(cat(1,min(trigger),min(force))) max(cat(1,max(trigger),max(force)))]);
grid on
filename=sprintf('force_data.eps');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 10 2]);
export_fig(filename)

%detect triggers
[tv,ti]=findpeaks(trigger,'MinPeakDistance',500);

%get force per cycle
force_sorted=zeros(size(ti,1)-1,sr*T);


for i=1:size(ti,1)-2        %don't include last cycle
    force_sorted(i,:)=force(ti(i):ti(i)+sr*T-1);
end

%get mean + low and up boundary (+/- std)
force_mean=mean(force_sorted,1);
force_std=std(force_sorted,0,1);
if min(force_mean)<0
    force_mean=force_mean+abs(min(force_mean));
elseif min(force_mean)>0
    force_mean=force_mean-abs(min(force_mean));
end

mvc_line=ones(size(force_mean))*MAX_force;

%create time array for plots
t=0:1/sr:T;
t(end)=[];


%find frame with peak force (according to numphases in mr acqusition)
[~,idx]=max(force_mean);
frame=round(idx/(sr*T)*numphases);

if permute_flag==1
% permuting to have max force in the middle of the cycle 
% permutation according to number of frames 
% not according to force sample rate
permute_ind=round(((sr*T)/numphases)*(numphases/2-frame));
force_mean=circshift(force_mean,permute_ind,2);
force_std=circshift(force_std,permute_ind,2);
end

force_high=force_mean+force_std;
force_low=force_mean-force_std;

%save force plot
h=figure;
set(gcf,'Visible','off');
plot(t,mvc_line,'--m','LineWidth', 1);
hold on
plot(t,force_mean,'-b','LineWidth', 2);
plot(t,force_high,'-r','LineWidth', 1);
plot(t,force_low,'-g','LineWidth', 1);

legend('MVC', 'Mean', 'Upper Bound', 'Lower Bound', 'Location','east')
ylabel('Voltage [V]');
xlabel('Time [s]');

axis([0 T 0-max(force_std) MAX_force*1.1])
axis square
grid on
filename=sprintf('force_average.eps');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.5 2.5 5 5]);
export_fig(filename)



%min force index
force_samplerate=size(force_mean,2)/numphases;
[~,point_idx]=min(force_mean(:));
force_ref=round(point_idx/force_samplerate);

%save data to *.mat
filename='force_data.mat';
save (filename,'MAX_force')
save (filename,'force_mean','-append')
save (filename,'force_std','-append')
save (filename,'frame','-append')
save (filename,'force_ref','-append')


end