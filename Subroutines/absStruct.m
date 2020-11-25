function [T]=absStruct(A)
fnames=fieldnames(A);

for i=1:size(A,2)
    for j=1:8
      
       A(i).(fnames{j})=abs(A(i).(fnames{j}));
    
    end
end

T=A;

