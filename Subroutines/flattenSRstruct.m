function a=flattenSRstruct(a)
% this function flattens SR array to make all fileds size 1

info = fieldnames(a);




    for i=1:size(info,1)
    
        t=numel(a(1).(info{i}));
    
        if t==3
            for l=1:3
                for j=1:size(a,2)
                a(j).([info{i},num2str(l)])=a(j).(info{i})(1,l);
                end
            end
            
            a=rmfield(a,info(i));
        elseif t==9
           
            for l=1:3
                for m=1:3
                    for j=1:size(a,2)
                     a(j).([info{i},num2str(l),num2str(m)])=a(j).(info{i})(l,m);
                    end
                end
            end

            a=rmfield(a,info(i));
        end
            



    end

end