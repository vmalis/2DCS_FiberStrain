function [index]=fiber_end(initial_x,initial_y,x,y)

%==========================================================================
% This function returns the index of a point which with respect to the 
% coordinates of initial point has the smallest angle 
% with the positive x axis
%
% INput:    x - array, n - number of the element
%
% OUTput:   value of the n-th largest element and it's index
%_________________________________________________________________________
%
% written by Vadim Malis
% 12/14 at UCSD RIL 
%==========================================================================

angle = zeros (size(x,1),1);

for i=1:size(x)
   
        v1 = [1;0];
        v2 = [x(i)-initial_x;initial_y-y(i)];

    angle(i) = atan2(abs(det([v1,v2])),dot(v1,v2))*180/pi;
    
end

[~,index]=min(angle);
