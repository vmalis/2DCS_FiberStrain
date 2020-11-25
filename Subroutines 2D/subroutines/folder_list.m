function [folder_list,number_of_folders] = folder_list(path)


% This function gives the list of folders excluding dummies (for MAC)
%   
%   [list of folders, number of folders]=folder_list(path to the folder)
%   
%_________________________________________________________________________
%
% written by Vadim Malis
% 12/14 at UCSD RIL 

cd(path)
path=dir;

for k = length(path):-1:1
    % remove non-folders
    if ~path(k).isdir
        path(k) = [ ];
        continue
    end

    % remove folders starting with .
    fname = path(k).name;
    if fname(1) == '.'
        path(k) = [ ];
    end
end

folder_list=path;
number_of_folders=length(path);
end