%% mean of a structure
function [M,SD] = mean_struct(b)


    fnames=fieldnames(b);

    for j=1:8% size(fnames,1)
    
        TEMP=[];
        for k=1:size(b,2)
            temp=getfield(b,{k},fnames{j});
            TEMP=cat(4,TEMP,temp);
        end
        
        mean_val=squeeze(mean(abs(TEMP),4));
        sd_val = std(abs(TEMP),0,4);
    
        M.(fnames{j})=mean_val;
        SD.(fnames{j})=sd_val;
    end
    
    
    for j=9:size(fnames,1)
    
        TEMP=[];
        for k=1:size(b,2)
            temp=getfield(b,{k},fnames{j});
            TEMP=cat(4,TEMP,temp);
        end
        
        mean_val=squeeze(mean((TEMP),4));
        sd_val = std(TEMP,0,4);
    
        M.(fnames{j})=mean_val;
        SD.(fnames{j})=sd_val;
    end



end