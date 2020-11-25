function []=strain_plot(data,stdev,units,suffix,limits,dt)

    
color=[0.9,0,0;0,0.7,0;0,0,0.95];
t = [0:dt:dt*(size(data,1)-1)];
t_int = [0:dt/10:dt*(size(data,1)-1)];

        figure
        hold on
        intData = interp1(1:size(data,1),data(:,1),1:0.1:size(data,1),'spline');
        p1=plot(t_int,intData,'-','Color',color(:,1));
        p12=plot(t,data(:,1),'o','Color',color(:,1));
        intData = interp1(1:size(data,1),data(:,2),1:0.1:size(data,1),'spline');
        p2=plot(t_int,intData,'-','Color',color(:,2));
        p22=plot(t,data(:,2),'o','Color',color(:,2));
        intData = interp1(1:size(data,1),data(:,3),1:0.1:size(data,1),'spline');
        p3=plot(t_int,intData,'-','Color',color(:,3));
        p32=plot(t,data(:,3),'o','Color',color(:,3));
        
        ylabel(units,'Interpreter','latex','FontSize',16)
        filename=[suffix,'.eps'];

        ylim(limits)
        


        errorbar(t,data(:,1),stdev(:,1),'o','Color',color(:,1));
        errorbar(t,data(:,2),stdev(:,2),'o','Color',color(:,2));
        errorbar(t,data(:,3),stdev(:,3),'o','Color',color(:,3));
      
        xlabel('$time \quad [\mathrm{s}]$','Interpreter','latex','FontSize',16)
        %title(ID,'Interpreter','latex');

        xlim([0 dt*(size(data,1)+1)])
        set(gcf,'color','w');
        set(gca,'TickLabelInterpreter','latex','FontSize',14);
        [~,objects] = legend([p12,p22,p32],{'$60 \% \, \mathrm{MVC}$','$40 \% \, \mathrm{MVC}$','$30 \% \, \mathrm{MVC}$'},'Interpreter','latex','FontSize',16,'location','best');
         objects(4).LineStyle = '-';
        objects(6).LineStyle = '-';
        objects(8).LineStyle = '-';
        
        export_fig(filename,'-eps');
        hold off
        close
        
        
end