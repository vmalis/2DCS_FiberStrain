%=========================================================================
%   Step0 for 2d Strain Rate Analysis Toolbox (master script)
%
%       part of 2d Strane Rate Toolkit
%=========================================================================
%
% This script
%     a) creates appropriate folder structure
%
% INput: user's input
%
% OUTput: created folder structure
%   
% written by Vadim Malis
% 01/15 at UCSD RIL
%==========================================================================


%% VEPC DATA

% Root Folder of study
window_info=sprintf('Choose Data_VEPC folder');
PathName = uigetdir('~/Desktop',window_info);
cd(PathName)
mkdir('DATA_VEPC');
cd('DATA_VEPC')

% Input subject initials (should be comma separateed, lower case will 
% be converted to the upper case)
x = inputdlg('Enter subject initials seperated with comma','Sample', [1 50]);
Subject_Initials=textscan(x{1},'%s','Delimiter',',');
Subject_Initials=upper(Subject_Initials{1});
f = figure('CloseRequestFcn',@ediTable);


% Column names and column format
columnname = {'Subject Initials','ISO_pre','ISO_post',...
    'ISO_post-rehab','ECC_pre','ECC_post','ECC_post-rehab'};
columnformat = {'char','logical','logical','logical','logical',...
    'logical','logical'};

% Define the data used in uitable; default values are true
s=num2cell(true(size(Subject_Initials,1),6));
d=cat(2,Subject_Initials,s);

% Create the uitable
t = uitable('Data', d,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', [false true true true true true true],...
            'RowName',[],...
            'ColumnWidth',{100,80,80,100,80,80,100});

% Set width height and position
ScreenSize=get(0,'Screensize');
tableextent = get(t,'Extent');
oldposition = get(t,'Position');
tnewposition = [oldposition(1) oldposition(2) tableextent(3) tableextent(4)];
set(t, 'Position', tnewposition);
set(f, 'Position',[ScreenSize(3)/2-(tableextent(3)+80)/2 ...
    ScreenSize(4)/2-(tableextent(4)+80)/2 ...
    tableextent(3)+40 tableextent(4)+40]);

% put tags to check selection wit ediTable function
set(f,'Tag', 'tableFigure')
set(t,'Tag','table')

waitfor(f)

% cell of selected studies, removing initials
selected_data(:,1)=[];
columnname(1)=[];

clearvars -except Subject_Initials selected_data columnname

% got selected study elements
[~,j]=ind2sub(size(selected_data),find([selected_data{:,:}]));


% Create study folders
studyfolder=ismember(1:size(selected_data,2),j);

for k=1:size(selected_data,2)

    if studyfolder(k)==1
        mkdir(columnname{k});
    end
end


% Create subject folders

for k=1:size(columnname,2)

subjectfolder=find([selected_data{:,k}]);
cd(columnname{k})

    for l=subjectfolder
        mkdir(Subject_Initials{l,1});
    end

cd ..
    
end


%% move each study modality to its own folder

fldr_list=folder_list(pwd);

mkdir('ISO')
mkdir('ECC')

for i=1:size(fldr_list,1)

    if fldr_list(i).name(1:3)=='ECC'
        movefile(fldr_list(i).name, [pwd '/ECC/']);
    
    elseif fldr_list(i).name(1:3)=='ISO'
    movefile(fldr_list(i).name, [pwd '/ISO/']);
                  
    else
    rmdir (fldr_list(i).name)        
    
    end
    
end


%% ------------------------------------------------------------------------
% Fiber DATA

% folders for fiber data (each dynamic study should have High-Res for fiber 
% coordinate)

cd ..
mkdir('DATA_FIBER');
cd('DATA_FIBER')

% Create subject folders
for k=1:size(Subject_Initials,1)

mkdir(Subject_Initials{k,1});    
cd(Subject_Initials{k,1})

    % per each study: just numbers, names doesn't matter, header ids are
    % used later
    
    for m=1:sum([selected_data{k,:}])
        mkdir(num2str(m))
    end
    
cd ..

end 