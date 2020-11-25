function [value,index]=max_n(x,n)

% This function returns the n-th largest element in the array
%
% INput:    x - array, n - number of the element
%
% OUTput:   value of the n-th largest element and it's index
%_________________________________________________________________________
%
% written by Vadim Malis
% 12/14 at UCSD RIL 

    if n>size(x(:));
    end

[x_sorted,idx]=sort(x,'descend');

value=x_sorted(n);
index=idx(n);
    