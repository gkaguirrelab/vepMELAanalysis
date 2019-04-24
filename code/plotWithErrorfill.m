function []=plotWithErrorfill(xdata,ydata,CI,edgecolor,fillcolor,markeredge,markerface)
    CI=cat(2,CI(1,:),fliplr(CI(2,:)));
    X=cat(2,xdata,fliplr(xdata));
    fill(X,CI,fillcolor,'EdgeColor',edgecolor);
    plot(xdata,ydata,'o','MarkerEdgeColor',markeredge,'MarkerFaceColor',markerface)
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
end