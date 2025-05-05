%% Example: fit Pearson's data with a straight line, using covariances given by York
% Author: anna.charvatovacampbell@cmi.gov.cz
% Created: 2025-04-04
 
clear
close all
addpath('../../')

% Load data
xydata =  dlmread("Pearson-data");
x = xydata(:,1);
y = xydata(:,2);

xcov = dlmread("York-xcov");
ycov = dlmread("York-ycov");

Uxx = diag(xcov);
Uyy = diag(ycov);

data = {x, y};
U      = {Uxx []; [] Uyy};

%alternatively
%load("PearsonYork.mat");

% Define function to be fitted as a constraint
fun = @(mu,beta) beta(1)*mu{1}+beta(2) - mu{2};
% Derivative of fun with respect to mu
options.funDiff_mu = @(mu, beta) {beta(1)*ones(size(mu{1})), -ones(size(mu{2}))};
% Derivative of fun with respect to beta
options.funDiff_beta = @(mu, beta) [mu{1}, ones(size(mu{1}))];


% Initial estimate of parameters and true values
beta0 = [1,1]';
mu0 = data;



% Define options
options.method = "oefpilrs2";
options.isPlot = false;
options.verbose = false;
options.maxit = 100;

% Fit the data
result = OEFPIL(data,U,fun,mu0,beta0,options);

% Check convergence
if result.iter == options.maxit
    fprintf("Warning: OEFPIL did not converge \n.");
end


% Print results
fprintf("Parameters y = beta_1*x + beta_2   \n");
fprintf("\t Best estimate \t Uncertainty\n");
for i=1:length(beta0)
    fprintf("beta_%d \t %g \t %g \n", i, result.beta(i), result.ubeta(i));
end

fprintf("\n");

fprintf("True values\n");
fprintf("\t Best estimate \t Uncertainty\n");
for i=1:length(result.mu(:,1))
    fprintf("mux_%d \t %g \t %g \n", i, result.mu(i,1), result.umu(i,1));
end
for i=1:length(result.mu(:,2))
    fprintf("muy_%d \t %g \t %g \n", i, result.mu(i,2), result.umu(i,2));
end

% Plot the results

x  = data{1};
y = data{2};
ux = sqrt(diag(U{1,1}));
uy = sqrt(diag(U{2,2}));

a = result.beta(1);
b = result.beta(2);

ua = result.ubeta(1);
ub = result.ubeta(2);
uab = result.Ubeta(1,2);


hold on;
xlabel("x");
ylabel("y");
title({"Straight line fit", "Pearson's data with York's weights"})

% calculate oefpil curve:
xx = linspace(min(x)-0.1*(max(x)-min(x)), max(x)+0.1*(max(x)-min(x)), 1000);
yy = fun({xx, zeros(size(xx))}, result.beta);
% calculate oefpil curve bands:
dyydbeta = options.funDiff_beta({xx',zeros(size(xx'))}, result.beta);
uyy = sqrt(diag(dyydbeta*result.Ubeta*dyydbeta'))';
% calculate polyfit function - least-square-error fit:
[pnls, snls] = polyfit(x,y,1);
[yynls, uyynls] = polyval(pnls, xx, snls);

% plot fills as first because fills should go to the bottom and plot rest over
% the fills
fill([xx, flip(xx)], [yynls-2*uyynls, flip(yynls+2*uyynls)], 'm', 'facealpha', 0.1, 'edgealpha', 0);
fill([xx, flip(xx)], [yy-2*uyy, flip(yy+2*uyy)], 'g', 'facealpha', 0.2, 'edgealpha', 0);
% plot data:
if isOctave
    errorbar(x, y, 2*ux, 2*uy, '~>.b');
else
    errorbar(x, y, 2*uy, 2*uy, 2*ux, 2*ux, '.b');
end
% plot oefpil curve and bands:
hoefpil = plot(xx, yy, 'g-', 'linewidth', 2);
                    % plot(xx, yy+2*uyy, 'g--', 'linewidth', 2);
                    % plot(xx, yy-2*uyy, 'g--', 'linewidth', 2);
% plot polyfit and bands:
hnls = plot(xx, yynls, 'm-');
                    % hnls = plot(xx, yynls - 2*uyynls, 'm--');
                    % hnls = plot(xx, yynls + 2*uyynls, 'm--');
% save plot to file:
legend([hnls, hoefpil], {'polyfit (least-square-error, )', 'OEFPIL'}, 'Location', 'northeast');
print('straight_line', '-dpng', '-r300')  % 300 DPI PNG
