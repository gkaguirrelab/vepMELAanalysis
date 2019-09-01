function []=plotWithXYErrorbars(xdata,ydata,xCI,yCI,markeredge,markerface)
    lower_errorX=xdata-xCI(1,:);
    upper_errorX=xCI(2,:)-xdata;
    lower_errorY=ydata-yCI(1,:);
    upper_errorY=yCI(2,:)-ydata;
    errorbar(xdata,ydata,lower_errorY,upper_errorY,lower_errorX,upper_errorX,'o','Color',markeredge,'MarkerEdgeColor',markeredge,'LineWidth',2,'MarkerFaceColor',markerface)
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
end