function [means,diffs,meanDiff,CR,linFit] = BlandAltman(var1, var2, flag,plot_title,X_lim,Y_lim,X_lable,Y_lable)
    
    
    if (nargin<1)
            display('Not enough input parameters!');
    end
    
    if nargin==2
        flag = 0;
    end
    
    means = mean([var1;var2]);
    diffs = var1-var2;
    
    meanDiff = mean(diffs);
    sdDiff = std(diffs);
    CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff]; %%95% confidence range
    
    linFit = polyfit(means,diffs,1); %%%work out the linear fit coefficients
    
    %%%plot results unless flag is 0
    if flag ~= 0
        plot(means,diffs,'ko','MarkerSize',10,'MarkerFaceColor','k')
        hold on
        if flag > 1
            plot(X_lim, ones(1,2).*CR(1),'r--','LineWidth',2); %%%plot the upper CR
            plot(X_lim, ones(1,2).*CR(2),'r--','LineWidth',2); %%%plot the lower CR
            plot(X_lim, zeros(1,2),'k'); %%%plot zero
            plot(X_lim, meanDiff*ones(1,2),'b','LineWidth',5); %%%plot zero
             
        end
        if flag > 2
            plot(means, means.*linFit(1)+linFit(2),'k--'); %%%plot the linear fit
        end
        
    end
    
    
    ax = gca;
    set(gcf,'color','w');
    set(gca,'TickLabelInterpreter','latex','FontSize',18);
    ylim(Y_lim);
    xlim(X_lim);
    xlabel(X_lable,'Interpreter','latex','FontSize',26)
    ylabel(Y_lable,'Interpreter','latex','FontSize',26)
    title(plot_title,'Interpreter','latex','FontSize',34)

end