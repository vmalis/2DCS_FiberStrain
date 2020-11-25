function c=meanANDsd(a,b)

    oldnames = fieldnames(b);

    for k=1:3
        for i=1:size(oldnames,1)
            newname=[oldnames{i},'_sd']; 
            oldname=oldnames{i};
            c(k).(oldname)=a(k).(oldnames{i});   
            c(k).(newname)=b(k).(oldnames{i});   
        end
    end    
end