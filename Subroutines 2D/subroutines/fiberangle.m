function [angle]=fiberangle(x1,y1,x2,y2)

%==========================================================================
%  Subroutine to calculate fiber angles using start- and endpoints info
%==========================================================================
%
% This subroutine is a part of 2D Strain rate Toolkit and is used by master
% script STEP2_ROIs. Subroutine calculates angles between fiber and
% positive x-axis
%
%
% INput:        x-start, y-start, x-end, y-end (column or row vector)
%
% OUTput:       angles (row vector [degrees])
%
%--------------------------------------------------------------------------
% written by Vadim Malis
% 12/14 at UCSD RIL
%==========================================================================


angle=zeros(size(x1));

    for i=1:size(x1)

        angle(i)=180/pi()*acos((x2(i)-x1(i))/sqrt((x2(i)-x1(i))^2+(y2(i)...
            -y1(i))^2));

    end
    
    
end
