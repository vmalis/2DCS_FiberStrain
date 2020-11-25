function B=stdfilter(A)

%=====================================================================
%subroutine takes an array of data, and outputs average value after
%removing entries that are more than 2 stdev away from the average

%part of ulls analysis

% written by Vadim Malis
% 05/16 at UCSD RIL
%=====================================================================

A_aver=nanmean(A,1);
index=[];

mmax=A_aver+nanstd(A,1);
mmin=A_aver-nanstd(A,1);

for i=1:size(A,1)
    
        if A(i)>mmax || A(i)<mmin
            index=[index,i];
        end

end

A(index)=[];
B=nanmean(A);
end