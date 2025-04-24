%% Example: fit Pearson's data with a straight line, using covariances given by York
% Author: anna.charvatovacampbell@cmi.gov.cz
% Created: 2025-04-04
 
clear
close all
addpath('../../')

% Load data
load("PearsonYork.mat");

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

%---------------------- PLOTTING -----------

x=data{1};
y=data{2};
ex = sqrt(diag(U{1,1}));
ey = sqrt(diag(U{2,2}));

a=result.beta(1);
b=result.beta(2);
xstar = result.mu(:,1);
ystar = result.mu(:,2);

ua = result.ubeta(1);
ub = result.ubeta(2);
uab= result.Ubeta(1,2);


hold on;
xlabel("x");
ylabel("y");
title("Straight line fit,\n Pearson's data with York's weights")
errorbar(x, y, 2*ey, "bo");

for i = 1:length(x)
        line([x(i) - 2*ex(i), x(i) + 2*ex(i)], [y(i), y(i)], ...
            'Color', 'blue', 'LineStyle', '-');
end


tt = linspace(min(x)-0.1*(max(x)-min(x)), max(x)+0.1*(max(x)-min(x)), 1000);

hoefpil=plot(tt, a*tt+b, "g-");
plot(tt, a*tt+b+2*sqrt(tt.^2*ua^2 + ub^2 + 2*tt*uab), "k--");
plot(tt, a*tt+b-2*sqrt(tt.^2*ua^2 + ub^2 + 2*tt*uab), "k--");

pnls=polyfit(x,y,1);
hnls=plot(tt,polyval(pnls,tt), "m");

legend([hoefpil, hnls], {'OEFPIL', 'NLS'}, 'Location', 'northeast');
print('straight_line', '-dpng', '-r300')  % 300 DPI PNG
