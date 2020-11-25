% function to create plots for the paper
function plots_4paper(data_raw,stdev,~,ylbl,filename)

data=data_raw;

color=[0.9,0,0;0,0.7,0;0,0,0.95];
lgnd={'Proximal','Middle','Distal'};

figure
hold on
intData = interp1(1:size(data,1),data(:,1),1:0.1:size(data,1),'spline');
p1=plot(1:0.1:size(data,1),intData,'-','Color',color(:,1));
p12=plot(data(:,1),'o','Color',color(:,1));
%eb1=errorbar(data(:,1),stdev(:,1),'o','Color',color(:,1));

intData = interp1(1:size(data,1),data(:,2),1:0.1:size(data,1),'spline');
p2=plot(1:0.1:size(data,1),intData,'-','Color',color(:,2));
p22=plot(data(:,2),'o','Color',color(:,2));
%eb2=errorbar(data(:,2),stdev(:,2),'o','Color',color(:,2));
        
intData = interp1(1:size(data,1),data(:,3),1:0.1:size(data,1),'spline');
p3=plot(1:0.1:size(data,1),intData,'-','Color',color(:,3));
p32=plot(data(:,3),'o','Color',color(:,3));
%eb3=errorbar(data(:,3),stdev(:,3),'o','Color',color(:,3));

legend([p1,p2,p3],lgnd)

xlabel('frame number')
ylabel(sprintf(ylbl));
%title(sprintf(ttl));
xlim([0 22.5])
set(gcf,'color','w');
export_fig(filename,'-eps');

hold off
end