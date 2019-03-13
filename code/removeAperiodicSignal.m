function []=removeAperiodicSignal(ttf,varargin)

%% Parse input
p = inputParser;
p.addParameter('TemporalFrequency',[1.625 3.25 7.5 15 30 60 90 120 180],@isnumeric);
p.addParameter('plot_all',false,@islogical);

p.parse(varargin{:});

ttf_M=ttf.ttf_M;
L=@(x0,xdata)x0(1)-log10(x0(2)+(xdata.^x0(3)));

%% Find (to later remove) power spectra at the fundamental and harmonic peaks of the stimuli
temp3=[];
for z=1:length(p.Results.TemporalFrequency)
    temp=find(ttf.f>TemporalFrequency(z));
    temp2=[temp(1)-1 temp(1)];
    temp3=cat(2,temp3,temp2);
end

%% Fit aperiodic curve
for x=1:size(ttf_M,1)
    xdata=ttf.f;
    xdata2=ttf.f;
    xdata2(:,temp3)=NaN;
    xdata2=rmmissing(xdata2);
    
    ydata=ttf_M(x,:);
    ydata2=ttf_M(x,:);
    ydata2(:,temp3)=NaN;
    ydata2=rmmissing(ydata2);
    
    x0=[log10(ydata(2)) 1 0.5];
    ub=[log10(ydata(2))+0.5 100 2];
    lb=[log10(ydata(2))-2 0.1 0.1];
    
    temp=find(xdata>200);
    L_fit(x,:)=lsqcurvefit(L,x0,xdata2(1:temp(1))',log10(ydata2(1:temp(1)))',lb,ub);
    LL=linspace(xdata(1),xdata(end),length(xdata));
    L_curve(x,:)=L(L_fit,LL);
    
    
    if p.Results.plot_all==1
        figure
        hold on
        plot(xdata,log10(ydata),'-k')
        plot(LL,L_curve(x,:),'-r')
        title(num2str(L_fit(x,:)))
        ax=gca;
        ax.TickDir='out';
        ax.Box='off';
        ax.XLim=[min(xdata) 100];
        hold off
        pause
    end
end

end
