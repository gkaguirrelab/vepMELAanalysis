function []=plotWithErrorbars(xdata,ydata,CI,markerline,markeredge,markerface)
    lower_error=ydata-CI(1,:);
    upper_error=CI(2,:)-ydata;
    errorbar(xdata,ydata,lower_error,upper_error,markerline,'MarkerEdgeColor',markeredge,'LineWidth',2,'MarkerFaceColor',markerface)
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
    ax.XLim=[min(xdata)-(0.9*min(xdata)) max(xdata)+(1.1*min(xdata))];
end