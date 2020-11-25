function ediTable(~,~)

%==========================================================================
%  Callback subroutine to save choice in the uitable
%==========================================================================
%   
%--------------------------------------------------------------------------
% written by Vadim Malis
% 12/14 at UCSD RIL
%==========================================================================


myfigure=findobj('Tag','tableFigure');
myData=get(findobj(myfigure,'Tag','table'),'Data');
assignin('base','selected_data',myData)
delete(myfigure)


end