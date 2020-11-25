function quiver_color(x,y,u,v,cont,units,flag)

%==========================================================================
%  Subroutine to plot vectorfield with colorcoded magnitude
%==========================================================================
%
% This subroutine is a part of 2D Strain rate Toolkit and is used by 
% supplementary script VF_plot. 
%
% INput:  [x,y]     meshgrid; [u,v] - vector components (not normalised)
%         cont      contour lines (definening color levels for colorcoding)
%         units     TeX input to put into colorbar description
%         flag      1 to have a colorbar, 0 to turn it off
%
%--------------------------------------------------------------------------
% written by Vadim Malis
% 02/15 at UCSD RIL
%==========================================================================

% Get current axis 
h=get(gcf,'CurrentAxes');

% Check sizes all agree in input data
if size(x)~=size(y); error('x and y sizes disagree'); end
if size(x)~=size(u); error('x and u sizes disagree'); end
if size(y)~=size(v); error('y and v sizes disagree'); end


% Remove masked grid points from the input by filling coordinates with NaN;
x(isnan(u))=NaN;
y(isnan(u))=NaN;
x=double(x);
y=double(y);
% Scale the vectors according to the reference arrow vector length based on
% the mean distance between grid points. This is a good measure, as it remains 
% constant for multiple plots using the same grid with different values.
% x1=abs(diff(x')); x2=abs(diff(x)); 
% y1=abs(diff(y')); y2=abs(diff(y));
% [~,z1] = cart2pol(x1,y1); [~,z2] = cart2pol(x2,y2);
% scalelength=min(mean(z1(~isnan(z1))),mean(z2(~isnan(z2))));
scalelength=1;

% Remove NaN values that will not be plotted
% and turn points into a row of coordinates
u=u(~isnan(x))';
v=v(~isnan(x))';
y=y(~isnan(x))';
x=x(~isnan(x))';

% Get the magnitude and direction of each vector
[th,z] = cart2pol(u,v);

% Check cont has at least one value
 if isempty(cont) || isempty(cont) || any(cont==0)
  error('cont must be non-zero and contain at least one value')
 end

 % Arrange contour values from min to max, and 
 % add an extra value for the purpose of processing
 cont=sort(abs(cont));
 cont(end+1)=cont(end)+1;

 % Set the colormap 
 hc=colormap(jet(length(cont)));

 for i=1:length(cont)

  if i==1
   mask=find(z<cont(i));
  elseif i==length(cont)
   mask=find(z>=cont(i-1));
  else
   mask=find(z<cont(i) & z>=cont(i-1));
  end

% Center vectors over grid points
  [u,v] = pol2cart(th(mask),scalelength);
  xstart=x(mask)-0.5*u;
  xend=x(mask)+0.5*u;
  ystart=y(mask)-0.5*v;
  yend=y(mask)+0.5*v;
  
%--------------------------------------------------------------------------
%Arrows can be added by uncommenting below
%  
% Set arrow size (1= full length of vector)
% arrow=0.33;  
  
% Get x coordinates of each vector plotted
%   lx = [xstart; ...
%        xstart+(1-arrow/3)*(xend-xstart); ...
%        xend-arrow*(u+arrow*v); ...
%        xend; ...
%        xend-arrow*(u-arrow*v); ...
%        xstart+(1-arrow/3)*(xend-xstart); ...
%        NaN(size(xstart))];
% 
% Get y coordinates of each vector plotted
%   ly = [ystart; ...
%        ystart+(1-arrow/3)*(yend-ystart); ...
%        yend-arrow*(v-arrow*u); ...
%        yend; ...
%        yend-arrow*(v+arrow*u); ...
%        ystart+(1-arrow/3)*(yend-ystart); ...
%        NaN(size(ystart))];
%--------------------------------------------------------------------------


lx = [xstart;xend];
ly = [ystart;yend];

% Plot the vectors
line(lx,ly,'Color',hc(i,:),'LineWidth',.01);

 end


% My fancy narow horizontal colorbar
% Use when you need a high quality images to be saved in file

if flag==1

h=colorbar('peer',gca,'location','SouthOutside');
caxis([cont(1),cont(end)])
set(h,'ylim',[cont(1),cont(end)]);
colorbar_scale=roundn((cont(end)-cont(1))/5,2)-cont(1);
colorbar_scale_r=roundn(colorbar_scale,1);
set(h,'YTick',[cont(1),cont(1)+colorbar_scale:colorbar_scale_r:cont(end)]);
xlabel(h,units,'FontSize', 8);
set(gca,'FontSize', 8);
cb_pos=get(h,'Position');
set(h,'Position',[cb_pos(1)*0.66,cb_pos(2)*.5,cb_pos(3)*1.32,cb_pos(4)*.3])
end

end

