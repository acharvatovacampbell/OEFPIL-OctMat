%% SThM minimalization

clear 
close all;

addpath('../../')


%% Set up data

caldata = dlmread("XY.dat");
x = caldata(2:end, 2);
ux = caldata(2:end, 3);
y = caldata(2:end, 4);
uy = caldata(2:end, 5);

unknown = dlmread("Y.dat");
z = unknown(2:end,2);
uz = unknown(2:end,3);

Ux = diag(ux.^2);
Uy = diag(uy.^2);
Uz = diag(uz.^2);

data = [x', y', z']';
udata = [ux', uy', uz']';

nx = length(x);
nz = length(z);

U = diag(udata.^2);


%% Set the function of parameter constraints

mu0 = data;
beta0=[[1,1,1], zeros(size(z'))]';
np =3;

fun = @(mu,beta) funSThM(mu,beta);
options.numDiffMethod = 'LynnesMoller';

%% Other options

options.tol = 1e-10;
options.method = 'oefpil';
options.criterion = 'parameterdifferences';
options.isEstimatedVariance=false;
options.maxit = 1000;
options.isPlot = false;

%% Run the OEFPIL algorithm

options.q = nx + nz;
result = OEFPIL(data,U,fun,mu0,beta0,options);

disp(['Number of iterations  ', num2str(result.iter)])
fprintf("\n");

disp('Parameters of calibration curve ');
fprintf(" \t estimate \t uncertainty\n");
abc=["a","b", "c"];
for i = 1:np
    fprintf('%s %g %g\n', abc(i), result.beta(i), result.ubeta(i))
end

fprintf("\n");
disp('Conductivity of unknown samples')
fprintf("sample \t estimate \t uncertainty\n");
unknownname=["A", "B", "D", "E"];
for i = 1:nz
    fprintf('%s %g %g\n', unknownname(i), result.beta(np+i), result.ubeta(np+i))
end

% Non-linear fitting if available
have_lsq = exist('lsqcurvefit');
if not(have_lsq)
    if isOctave % octave environment
        try
            pkg load optim;
            have_lsq = true;
        catch errmsg1
            disp('Non-linear fitting in this example requires GNU Octave package `optim`. You may consider to install it from Octave-forge or proceed without it.')
        end % try
    else % matlab environment
        error('Non-linear fitting in this example requires the function `lsqcurvefit` from Matlab package `Optimization Toolbox`. Either proceed without it, or use open-source GNU Octave.');
    end %if
end 

if have_lsq
   p0 = [1,1,1];
   if isOctave
        ff = @(x, p) p(1) * x./(p(2) + x) + p(3);
        [~, pnls, ~, ~, ~, covpar, ~, ~, ~, ~] = leasqr(x,y, p0, ff);
        upnls = sqrt(diag(covpar)); 
   else
        ff = @(p, x) p(1) * x ./(p(2) + x) + p(3);
        [pnls, resnorm, residuals, exitflag, output, lambda, jacobian]= lsqcurvefit(ff, p0, x, y);
        dof = length(y) - np;
        mse = resnorm / dof;
        covpar = mse * inv(jacobian' * jacobian);
        upnls = sqrt(diag(covpar)); 
   endif
   disp("\n\nNon-linear least squares");
   disp('Parameters of calibration curve ');
   for i = 1:np
        fprintf('%s %g %g \n', abc(i), pnls(i), upnls(i))
   end

   fprintf("\n");
   invval = pnls(2)*(pnls(3)-z)./(z-pnls(1)-pnls(3));
   
   disp('Conductivity of unknown samples')
%   dydp = [ pnls(2)*(pnls(3)-z)./(z-pnls(1)-pnls(3)).^2, (pnls(3) - z)./(z-pnls(1)-pnls(3)), 
%            pnls(3)./(z-pnls(1)-pnls(3)) + pnls(2)*(pnls(3)-z)./(z-pnls(1)-pnls(3)).^2 ];
   dydp1 = pnls(2)*(pnls(3)-z)./(z-pnls(1)-pnls(3)).^2;
   dydp2 = (pnls(3) - z)./(z-pnls(1)-pnls(3));
   dydp3 = pnls(3)./(z-pnls(1)-pnls(3)) + pnls(2)*(pnls(3)-z)./(z-pnls(1)-pnls(3)).^2;
   dydp = [ dydp1, dydp2, dydp3];

   dydz = -pnls(2)/(z-pnls(1)-pnls(3)) - pnls(2)*(pnls(3)-z)./(z-pnls(1)-pnls(3)).^2;

   uinvval2 = dydz^2*uz.^2 + diag(dydp*covpar*dydp');
   uinvval = sqrt(uinvval2);
   fprintf("sample \t estimate \t uncertainty\n");
   for i = 1:nz
     fprintf('%s %g %g \n', unknownname(i), invval(i), uinvval(i) )
   end
end

% Plot the results

a = result.beta(1);
b = result.beta(2);
c = result.beta(3);

errorbar(x, y, 2*ux, 2*uy, "~>.b");
hold on;

lt = log10(min(x))-0.04*(log10(max(x))-log10(min(x)));
ut = log10(max(x))+0.04*(log10(max(x))-log10(min(x)));
tt = linspace(lt, ut, 501)';
xx = 10.^(tt);
yy = result.beta(1)*xx./(result.beta(2)+xx)+ result.beta(3);
xlim([10^lt, 10^ut]);


hoefpil = plot(xx, yy, "g-");

dy=[ xx./(result.beta(2) + xx), -result.beta(1).*xx./(result.beta(2)+xx).^2, ones(size(xx))];
uy = diag(sqrt(dy*result.Ubeta(1:np,1:np)*dy'));

plot(xx, yy + 2*uy, "k--");
plot(xx, yy - 2*uy, "k--");

if have_lsq
   if isOctave
        hnls = plot(xx, ff(xx,pnls), "m");
   else
       hnls = plot(xx,ff(pnls, xx), "m");
   end
   legend([hoefpil, hnls], {'OEFPIL', 'NLS'}, 'Location', 'southeast');
else
   legend(hoefpil, 'OEFPIL', 'Location', 'southeast');
end

set(gca, 'XScale', 'log')
xlabel("thermal conductivity (W m^{-1} K^{-1})");
ylabel("Y measurand");
title("SThM calibration curve,\n ax/(b+x) +c ")
print('SThMcalibration', '-dpng', '-r300')  % 300 DPI PNG

