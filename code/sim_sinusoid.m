function [yF]=sim_sinusoid(t,freq,y)


% The sine-forming vector for one cycle in the measurement window
s = t./max(t) * 2 * pi;

% The Fourier basis. Note that the dimensions of X are tx2
X=[];
X(:,1) = sin(s * freq);
X(:,2) = cos(s * freq);

% The regression
b=regress(y,X);

% The modeled response
yF = X*b;

end