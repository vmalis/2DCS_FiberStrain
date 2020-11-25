function new_struct=struct_replace(old_stuct,fieldname,what_to_replace,with_what_to_replace)
%==========================================================================
%  Subroutine to replace find and replace data of string type in the struct
%==========================================================================
%
% This subroutine is a part of 2D Strain rate Toolkit and is used by 
% Step3_analysis. 
%
% INput:  struct
%         field to replace
%         2 arrays (1st - what to replace;  2nd - with what to replace)
%--------------------------------------------------------------------------
% written by Vadim Malis
% 02/15 at UCSD RIL
%==========================================================================

n=size(what_to_replace,2);
X={old_stuct.(fieldname)};

for i=1:n
index= find(strcmp(X, what_to_replace(i)));
[old_stuct(index).(fieldname)]=deal(with_what_to_replace(i));

end

new_struct=old_stuct;