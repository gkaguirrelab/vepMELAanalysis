% Generate primary settings as a function of eccentricity

% Define the eccentricity support
x = 1:1:30;

% Define a power law function with a y-offset
powerLaw = @(a,n,x) a+x.^n;

% Define the post-receptoral directions to examine
directions = {'S','LMinusM'};
receptorContrastIdx = [3,1,1];


for dd = 1:length(directions)
    
    % Obtain the modulation primaries for positive S as a function of
    % eccentricity
    for ee = x
        resultSet = designNominalSPDs('observerAgeInYears',32,'primaryHeadRoom',0.1,'fieldSizeDegrees',ee,'saveDir','');
        contrast(ee) = resultSet.(directions{dd}).positiveReceptorContrast(receptorContrastIdx(dd));
        RGB(ee,:)=resultSet.(directions{dd}).modulationPrimary;
    end
    
    % Obtain the a and n params for each of the primaries
    for pp = 1:3
        fitResult = fit(x',RGB(:,pp),powerLaw,'StartPoint',[0 0]);
        a(pp) = fitResult.a;
        n(pp) = fitResult.n;
    end
    
    % Plot the functions and fits
    plotLines = {'-r','-g','-b'};
    fitLines = {'*r','*g','*b'};
    figure
    hold on
    for pp = 1:3
        plot(x,RGB(:,pp),plotLines{pp})
        plot(x,powerLaw(a(pp),n(pp),x),fitLines{pp})
    end
    xlabel('Eccentricity [deg]')
    ylabel('Primary setting [0-1]')
    title(directions{dd});
    hold off
    
    % Report the functions
    fprintf(['Direction: ' directions{dd} '\n']);
    primaryLabels = {'R','G','B'};
    for pp = 1:3
        outline = [primaryLabels{pp} 'primary (deg) = ' num2str(a(pp)) ' + deg .^ ' num2str(n(pp)) '\n'];
        fprintf(outline);
    end
    fprintf('\n');
end

