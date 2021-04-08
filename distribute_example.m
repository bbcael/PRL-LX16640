clc; clear; close all;

E = 1e-8; % bulk epsilon
xi = 18.4152; % parameters
omega = 2.8145;
alpha =  1.9319;
% use the provided log-skew-normal fitting code to estimate these for a 
% given dataset; global values above. note xi is positive because we fit to 
% the difference from the smallest value (7e-14); this is only to get the 
% Owen's T-function code to work and has no bearing on the results. note 
% the fitting is also using the natural logarithm.

epsmin = log(7e-14); % lowest value considered (taken from lowest value
% in global dataset)
epsmax = log(1e-6); % highest value considered
eps = 0:.001:(epsmax-epsmin);

% generate cdf of log-skew-normal
lsncdf = normcdf((10.*log(eps)-xi)/omega)-2.*tfn((10.*log(eps)-xi)/omega,alpha);
% differentiate to get pdf
nspdf = diff(lsncdf);
% corresponding actual epsilon values
vals = exp(epsmin + eps);
vals = vals(2:end);
% total epsilon
tot = sum(vals.*nspdf);
% correct by E to get correct total
nrm = E./tot;
nspdf = nrm.*nspdf;
cont = nspdf.*vals;
clear eps lsncdf nrm tot nspdf;


semilogx(vals,cont,'k','linewidth',3) % cont(i) is now the amount of dissipation occurring, for each epilson 
% value in vals, where epsilon = vals(i)
axis([-Inf Inf -Inf Inf])
box on
set(gca,'fontsize',16,'ticklabelinterpreter','latex')
xlabel('$\varepsilon$','interpreter','latex')
ylabel('Contribution to bulk $\mathcal{E}$','interpreter','latex')
% note the magnitude of the y-value in the plot is resolution-dependent
% because sum(cont) = E
