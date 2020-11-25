%% devide structure by structure

function [S] = dev_struct(a,b)


    fnames=fieldnames(a);

    for j=1:1:size(fnames,1)
    
        
        S.(fnames{j})=a.(fnames{j})./b.(fnames{j});
        

    end



end