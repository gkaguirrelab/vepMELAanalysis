function []=plotWithErrorfill(xdata,ydata,CI,edgecolor,fillcolor,markeredge,markerface)
    CI=cat(2,CI(1,:),fliplr(CI(2,:)));
    X=cat(2,xdata,fliplr(xdata));
    fill(X,CI,fillcolor,'EdgeColor',edgecolor);
    if markeredge~='none'
        plot(xdata,ydata,'o','MarkerEdgeColor',markeredge,'MarkerFaceColor',markerface,'MarkerSize',10)
    end
    ax=gca;
    ax.TickDir='out';
    ax.Box='off';
end