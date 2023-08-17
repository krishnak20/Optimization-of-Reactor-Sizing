CA0 = [2 5 6 6 11 14 16 24];
Ca_in = 10;
Ca_out = 1;
CAf = [0.5 3 1 2 6 10 8 4];
Rate = [20 0.5 10 2 0.8 5 2.5 0.2];
T = [30 1 50 8 4 20 20 4];
XA = (CA0-CAf)./CA0;
V0 = 0.1;
rA = (CA0-CAf)./T;
Rate = 1./rA;
% Fit a spline curve to the data
pp = spline(CAf, Rate); % piecewise polynomial.

% Evaluate the spline curve at a fine grid of points
xgrid = linspace(min(CAf), max(CAf), 1000);
yfit = ppval(pp, xgrid);
% Plot the data points and the spline curve

figure;
scatter(CAf, Rate, 'filled');
hold on;
plot(xgrid, yfit, '-r', 'LineWidth', 2);
grid on;
xlabel('CAf (in mmol/m^3)');
ylabel('-1/rA (in min.m^3/mmol)');
title('Spline Interpolation of -1/rA vs CA');
legend('Data', 'Spline Curve');
f = @(x) ppval(pp, x);

% 1) For a single PFR
T_PFR = integral(f, Ca_out,Ca_in); 
Volume_PFR = V0*T_PFR;

% 2) For a single CSTR
 T_CSTR = (Ca_in-Ca_out)*(f(Ca_out));
 Volume_CSTR = V0*T_CSTR;

% 3) Two stirred tanks of any size
Ca_optimum = NaN;
rateVal = NaN;
for i = Ca_out:0.00001:Ca_in
    m = (f(Ca_out)-f(i))/(i-Ca_in);
    h = 1e-5; % a small number
    fx = ppval(pp, i);
    fxh = ppval(pp, i + h);
    dydx = (fxh - fx) / h;
    if abs((m - dydx)) < 1e-4  % Set tolerance for matching slopes
        Ca_optimum = i;
        rateVal = f(i);
        break;
    end
end

V1 = V0*(Ca_optimum-Ca_out)*f(Ca_out);
V2 = V0*(Ca_in-Ca_optimum)*f(Ca_optimum);
Volume_TwoStirred = V1+V2;

% 4) A combination of a PFR and a MFR

CaMin = fminbnd(f, Ca_out,Ca_in);
RateMin = ppval(pp, CaMin);
Vpfr = (integral(f, Ca_out,CaMin))*V0; 
Vcstr = (Ca_in-CaMin)*RateMin*V0; 
Volume_Combination = Vpfr+Vcstr;

% 5) A PFR with recycle
XAf = (Ca_in-Ca_out)/Ca_in;

for XAi = 0:0.001:XAf
    lowlimit = Ca_in*(1-XAf);
    upLimit = Ca_in*(1-XAi);
    LHS =  f(Ca_in*(1-XAi))*(XAf-XAi)*Ca_in;
    RHS = integral(f,lowlimit,upLimit);
    if abs(LHS - RHS) < 0.1 % set a tolerance for LHS - RHS
        conversion = XAi;
        break % exit the loop if a solution is found
    end
end

Ca_opt = Ca_in*(1-conversion);
RecyleRatio = (Ca_in-Ca_opt)/(Ca_opt-Ca_out);
Volume_recycle = (Ca_in-Ca_out)*f(Ca_opt)*V0;
VolumetricFlowRecycle = V0*RecyleRatio;


disp('All the Volumes are in m^3 and Concentrations are in mmol/m^3 \n')

fprintf('Volume of Single PFR is %f\n', Volume_PFR);
fprintf('Volume of Single CSTR is %f\n', Volume_CSTR);
fprintf('Volume of two stirred tank is %f\n', Volume_TwoStirred);
fprintf('Optimum Concentration is %f\n', CaMin);
fprintf('Volume of combination of a PFR and a MFR is %f\n', Volume_Combination);
fprintf('Volume of a PFR with recycle is %f\n', Volume_recycle);
fprintf('Optimum Recycle ratio of PFR is %f\n', RecyleRatio);
fprintf('Inlet Concentration in case of Recycle is %f\n', Ca_opt);



%CSTR
figure;
x1 = linspace(Ca_out, Ca_in, 1000);
y1 = ppval(pp, x1);
fill([Ca_out, Ca_out, Ca_in, Ca_in], [0, f(Ca_out), f(Ca_out), 0], 'blue', 'FaceAlpha', 0.3);
hold on;
plot(xgrid, yfit, '-k', 'LineWidth', 2);
grid on;
xlabel('CAf (in mmol/m^3)');
ylabel('-1/rA (in min.m^3/mmol)');
title('Shaded Area gives Residence Time for CSTR');
legend('Area', 'Spline Curve');

%PFR
figure;
x2 = linspace(Ca_out, Ca_in, 1000);
y2 = ppval(pp, x2);
area(x2, y2, 'FaceColor', 'blue', 'FaceAlpha', 0.3);
% line(x2, y2, 'Color', 'r', 'LineWidth', 2);
hold on;
plot(xgrid, yfit, '-k', 'LineWidth', 2);
grid on;
xlabel('CAf (in mmol/m^3)');
ylabel('-1/rA (in min.m^3/mmol)');
title('Shaded Area gives Residence Time for PFR');
legend('Area', 'Spline Curve');

% PFR and CSTR
figure;
x3 = linspace(Ca_out, CaMin, 1000);
y3 = ppval(pp, x3);
hold on;
area(x3, y3, 'FaceColor', 'blue', 'FaceAlpha', 0.3);
fill([CaMin, CaMin, Ca_in, Ca_in], [0, f(CaMin), f(CaMin), 0], 'red', 'FaceAlpha', 0.3);
% line(x3, y3, 'Color', 'r', 'LineWidth', 2);
plot(xgrid, yfit, '-k', 'LineWidth', 2);
xticks(CaMin);
grid on;
xlabel('CAf (in mmol/m^3)');
ylabel('-1/rA (in min.m^3/mmol)');
title('Left Shaded Area gives Residence Time of PFR, Right one gives CSTR');
legend('For PFR', 'For CSTR', 'Spline Curve');

% Two Stirred Tank

figure;
x4 = linspace(Ca_out, Ca_in, 1000);
y4 = ppval(pp, x4);
fill([Ca_out, Ca_out, Ca_optimum, Ca_optimum], [0, f(Ca_out), f(Ca_out), 0], 'blue', 'FaceAlpha', 0.3);
hold on;
fill([Ca_optimum, Ca_optimum, Ca_in, Ca_in], [0, f(Ca_optimum), f(Ca_optimum), 0], 'green', 'FaceAlpha', 0.3);
% line(x4, y4, 'Color', 'r', 'LineWidth', 2);
plot(xgrid, yfit, '-k', 'LineWidth', 2);
xticks(Ca_optimum);
grid on;
xlabel('CAf (in mmol/m^3)');
ylabel('-1/rA (in min.m^3/mmol)');
title('Left Shaded Area gives Residence Time of 1st CSTR, Right one gives 2nd CSTR');
legend('1st CSTR','2nd CSTR', 'Spline Curve');


%PFR with recycle
figure;
x5 = linspace(Ca_out, Ca_opt, 1000);
y5 = ppval(pp, x5);
area(x5, y5, 'FaceColor', 'green', 'FaceAlpha', 0.3);
% line(x5, y5, 'Color', 'r', 'LineWidth', 2);
% line(x5, zeros(size(y5)), 'Color', 'r', 'LineWidth', 2);
hold on;
plot(xgrid, yfit, '-k', 'LineWidth', 2);
xticks(Ca_opt);
grid on;
xlabel('CAf (in mmol/m^3)');
ylabel('-1/rA (in min.m^3/mmol)');
title('Shaded Area gives Residence Time for PFR with Recycle');
legend('Area', 'Spline Curve');
