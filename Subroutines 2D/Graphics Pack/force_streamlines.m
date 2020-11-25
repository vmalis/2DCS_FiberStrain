
t=[0:3/600:3];
t(end)=[];
min_force=min(force_mean);
force_mean=force_mean-min_force;

for i=1:22
hfig=figure('Color','black','Position', [100, 100, 900, 200])

x=round(600/22*i);

line(t(1:x),force_mean(1:x),'Color','g','LineWidth',2)
ax1 = gca;
ax1.XColor = 'white';
ax1.YColor = 'white';

ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation',...
           'right','Color','none');
set(ax2,'YTickLabel','')
set(ax2,'ytick',[])
set(ax2,'XTick',[1:22])
set(ax2,'xlim',[0,22]);
set(ax1,'ylim',[0,350]);
set(ax1,'xlim',[0,3]);
ax2.XColor = 'white';
ax2.YColor = 'white';

ax1.Color = 'black';

xlabel(ax1,'time [s]')
xlabel(ax2,'frame')
ylabel(ax1,'Force [N]')
line([x,x],[0,250],'Color','y');


filename=sprintf('force%2d', i);
export_fig(filename,'-jpg','-m4')
close
end



% 
% 
% hfig=figure('Color','white');
% 
% hl1=line(t,F_Y,'Color','b','LineWidth',2);
% hold on
% hl2=line(t,F_O,'Color','r','LineWidth',2);
% legend('Young','Senior');
% ax1 = gca;
% ax1_pos = ax1.Position;
% ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation',...
%            'right','Color','none');
% set(ax2,'YTickLabel','')
% set(ax2,'ytick',[])
% set(ax2,'XTick',[1:22])
% set(ax2,'xlim',[0,22]);
% set(ax1,'ylim',[0,400]);
% set(ax1,'xlim',[0,3]);
% 
% xlabel(ax1,'time [s]')
% xlabel(ax2,'frame')
% ylabel(ax1,'Force [N]')
% 
% 
% filename=sprintf('force');
% export_fig(filename,'-png','-m16')
% close
