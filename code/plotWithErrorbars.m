function []=plotWithErrorbars(xdata,ydata,CI,markerline,markeredge,markerface)
    lower_error=ydata-CI(1,:);
    upper_error=CI(2,:)-ydata;
    errorbar(xdata,ydata,lower_error,upper_error,markerline,'MarkerEdgeColor',markeredge,'LineWidth',1,'MarkerFaceColor',markerface,'MarkerSize',8)
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
end