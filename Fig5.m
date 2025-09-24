% Fig 5: multi-deme E and max Rj with empirical data from Veneto, Italy
clearvars; clc; close all; tic;

% Assumptions and notes
% - accounts for different starting times
% - optimal designs achieved by sampling from qR distributions
% - includes an E and max Rj metrics
% - no interactiion among demes and same serial intervals used

% Directory and where saving (or loading)
thisDir = cd; saveFol = 'data/veneto'; 
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Load incidence data and serial intervals
cd(saveFol);
% Read and count files with incidence data
files = dir('province*'); nDeme = length(files);
% Extract dates and new cases
Idata = cell(1, nDeme);
for i = 1:nDeme
    Idata{i} = readtable(files(i).name);
    % Initial data has nan so remove
    Idata{i} = Idata{i}(2:end, :);
end
cd(thisDir);

%% Format empirical data and obtain total infectiousness

% Stop time and full length of time
tdate = Idata{1}.date; tstop = find(tdate == '15/12/2020'); 
tlen = length(Idata{1}.date); nday = min(tstop, tlen);

% Extract dates and incidence
tday = 1:nday; tdate = tdate(1:nday); Ideme = zeros(nDeme, nday);
for i = 1:nDeme
    % Assume all files of same length nday
    Ideme(i, :) = Idata{i}.new_cases(tday);
    % Remove NaNs and replace with 0s
    Ideme(i, isnan(Ideme(i,:))) = 0;
    % Add smoothing (trailing)
    Ideme(i, :) = round(movmean(Ideme(i, :), [6 0]));
end

% Truncate time series so starts from first non-zero term
ndays = zeros(1, nDeme); istarts = zeros(1, nDeme);
for i = 1:nDeme
    istarts(i) = find(Ideme(i, :) > 0, 1, 'first');
    ndays(i) = nday - istarts(i) + 1;
end

% Define a standard serial interval distribution
wmean = 4.7; wvar = 2.9^2;
% Compose as a gamma distribution
scalePm = wvar/wmean; shapePm = wmean/scalePm;
wch = gamcdf(tday, shapePm, scalePm) - gamcdf(tday-1, shapePm, scalePm);

% Total infectiousness of each deme
Ldeme = Ideme; 
for j = 1:nDeme
    for i = 2:nday
        Ldeme(j, i) = sum(Ideme(j, i-1:-1:1).*wch(1:i-1));
    end
end

% Aggregate the epidemics and also Lam as serial interval fixed
Itot = sum(Ideme, 1); Ltot = sum(Ldeme, 1);

%% EpiFilter estimates of R, D and E numbers

% Grid limits and noise level
Rmin = 0.1; Rmax = 20; eta = 0.1;
% Uniform prior over grid of size m
m = 2000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Estimates variables for EpiFilter
Rm = zeros(nDeme, nday); Rl = Rm; Rh = Rm; qR = cell(1, nDeme);

% Smoothed estimates and distributions
for i = 1:nDeme
    % EpiFilter estimate from each deme
    [~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, ndays(i), p0,...
        Ldeme(i, istarts(i):end), Ideme(i, istarts(i):end));
    [~, Rl(i, istarts(i):end), Rh(i, istarts(i):end), Rm(i, istarts(i):end),...
        qR{i}] = runEpiSmoother(Rgrid, m, ndays(i), pR, pRup, pstate);
    
    disp(['Completed deme ' num2str(i) ' of ' num2str(nDeme)]);
end

% Aggregrate estimate over demes
[~, ~, ~, ~, pRL, pRupL, pstateL] = runEpiFilter(Rgrid, m, eta, nday, p0, Ltot, Itot);
[~, RLaml, RLamh, RLam, qRLam] = runEpiSmoother(Rgrid, m, nday, pRL, pRupL, pstateL);
clearvars('pstate', 'pstateL', 'pRL', 'pRupL', 'pR', 'pRup');

% Basic D and E optimal design means (no CIs)
Re_D = mean(Rm); Re_E = mean(Rm.^2)./Re_D;
% Prob > 1 for metrics
p1E = zeros(1, nday); p1D = p1E; p1R = p1D; p1Rjmax = p1D;

% D and E optimal distribution by sampling and max Rj
nsamps = 50000; RmD = zeros(1, nday); RlD = RmD; RhD = RmD; RmE = RmD; 
RlE = RmD; RhE = RmD; Rhmaxj = RmD; Rlmaxj = RmD; Rmmaxj = RmD;

for i = 1:nday
    % Individual distribution samples
    xDeme = zeros(nDeme, nsamps);
    for j = 1:nDeme
        if i >= istarts(j)
            xDeme(j, :) = datasample(Rgrid, nsamps, 'Weights', qR{j}(i-istarts(j)+1, :));
        end
    end
    % D and E optimal samples for this day
    Dsamp = mean(xDeme); 
    Esamp = mean(xDeme.^2)./Dsamp;
    % Sample maximum across demes
    Rjmaxsamp = max(xDeme);

    % Statistics of D designs
    RmD(i) = mean(Dsamp); 
    Dquants = quantile(Dsamp, [0.025, 0.975]);
    RlD(i) = Dquants(1); RhD(i) = Dquants(2);

     % Statistics of E designs
    RmE(i) = mean(Esamp);
    Equants = quantile(Esamp, [0.025, 0.975]);
    RlE(i) = Equants(1); RhE(i) = Equants(2);

    % Statistics of max of Rj metric
    Rmmaxj(i) = mean(Rjmaxsamp);
    Rjmaxquants = quantile(Rjmaxsamp, [0.025, 0.975]);
    Rlmaxj(i) = Rjmaxquants(1); Rhmaxj(i) = Rjmaxquants(2);

    % Prob of D > 1, E > 1, max Rj > 1
    p1D(i) = length(Dsamp(Dsamp > 1))/nsamps;
    p1E(i) = length(Esamp(Esamp > 1))/nsamps;
    p1Rjmax(i) = length(Rjmaxsamp(Rjmaxsamp > 1))/nsamps;
end

id1 = find(Rgrid <= 1, 1, 'last');
for i = 1:nday
    p1R(i) = 1 - sum(qRLam(i, 1:id1), 2)';
end


%% Figures on overall and individual deme performance

% Time to narrow down simulations on
tstart = find(tdate == '15/08/2020');
trange = [tstart tstop];

% Peak sizes of epidemics
figure('Position', [10 10 800 1000]);
subplot(4, 1, 1);
plot(tdate', Ideme', 'LineWidth', 2);
hold on;
plot(tdate', Itot', 'k', 'LineWidth', 2);
grid off; box off; hold off;
h = gca; h.YScale = 'log';
ylabel('$I_j(t)$', 'FontSize', fnt);
xlim(tdate(trange));

% Get ticks for middle panel
xt = h.XTick; xtlab = h.XTickLabel;
% Find time points of ticks 
tt = zeros(size(xt));
for i = 1:length(xt)-1
    tt(i) = find(tdate == xt(i));
end
tt(length(xt)) = tday(end)+1;

% Reproduction number estimates
subplot(4, 1, [2 3]);
hold on;
plotCIRaw(tday', Rmmaxj', Rlmaxj', Rhmaxj', 'c');
plotCIRaw(tday', RLam', RLaml', RLamh', 'b');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$\max R_j$', '', '$R$', '', '$E$');
legend('boxoff');
ylabel('$X(t)$', 'FontSize', fnt);
xlim(tday(trange)); ylim([0.7 2.1]);

% Probability of resurgence
subplot(4, 1, 4);
hold on; 
plot(tdate, p1Rjmax, '-', 'color', 'c', 'linewidth', 2);
plot(tdate, p1R, '-', 'color', 'b', 'LineWidth', 2);
%plot(tdate, p1D, '-', 'color', 'g', 'linewidth', 2);
plot(tdate, p1E, '-', 'color', 'r', 'LineWidth', 2);
plot(tdate, 0.5*ones(1, nday), '--', 'color', 'k', 'LineWidth', 1);
h = gca; h.YColor = h.XColor;
grid off; box off; hold off;
xlim(tdate(trange)); ylim([0 1.05]);
ylabel('P($X(t) > 1$)', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt);

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Individual deme plots
ax1 = zeros(1, nDeme); ax2 = ax1; grey3 = [0.9 0.9 0.9];
% Comparison of each deme estimates
figure('Position', [10 10 700 1100]);
for i = 1:nDeme
    ax1(i) = subplot(nDeme+2, 1, i);
    yyaxis('right');
    plot(tday, Ideme(i, :), 'Color',  'k', 'LineWidth', 2);
    h = gca; h.YColor = h.XColor;
    yyaxis('left');
    plotCIRaw(tday', Rm(i, :)', Rl(i, :)', Rh(i, :)', 'b');
    hold on; h = gca; h.YColor = h.XColor;
    plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
    grid off; box off; hold off; xlim([tstart tstop]); 
    ylabel(['$R_{' num2str(i) '}(t)$'], 'FontSize', fnt);
end

% Total and consensus R estimates
ax1(nDeme+1) = subplot(nDeme+2, 1, [nDeme+1 nDeme+2]);
yyaxis('right');
stairs(tday, Itot, 'Color', 'k', 'LineWidth', 2);
h = gca; h.YColor = h.XColor;
yyaxis('left');
plotCIRaw(tday', Rmmaxj', Rlmaxj', Rhmaxj', 'c');
hold on; h = gca; h.YColor = h.XColor;
plotCIRaw(tday', RLam', RLaml', RLamh', 'b');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
ylabel('$X(t)$', 'FontSize', fnt);
legend('', '$\max R_j$', '', '$R$','', '$E$');
legend boxoff;
xlabel('$t$ (days)', 'FontSize', fnt);
xlim([160 290]); linkaxes(ax1, 'y'); 

%% Publishable figure 5

% Middle cut point
tmiddle = find(tdate == '15/10/2020');
trange1 = [tstart tmiddle]; t1 = tstart:tmiddle;
trange2 = [tmiddle+1 tstop]; t2 = tmiddle+1:tstop;

figure('Position', [10 10 1000 800]);
% Epidemic trajectories
subplot(4, 2, 1);
plot(tdate(t1)', Ideme(:, t1)', 'g-', 'LineWidth', 2);
hold on; yymax = max(max(Ideme(:, t1)));
plot(tdate(t1)', yymax*Itot(t1)'/max(Itot(t1)), 'k-', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$I_j(t)$', 'FontSize', fnt);
xlim(tdate(trange1)); h = gca;
h.XTick = [tdate(t1(1)), tdate(ceil(median(t1))), tdate(t1(end))];

% Get ticks for middle panel
xt = h.XTick; xtlab = h.XTickLabel;
% Find time points of ticks 
tt = zeros(size(xt));
for i = 1:length(xt)-1
    tt(i) = find(tdate == xt(i));
end
tt(length(xt)) = tday(t1(end));

% Reproduction number estimates
subplot(4, 2, [3 5]);
hold on;
plotCIRaw(tday(t1)', Rmmaxj(t1)', Rlmaxj(t1)', Rhmaxj(t1)', 'c');
plotCIRaw(tday(t1)', RLam(t1)', RLaml(t1)', RLamh(t1)', 'b');
plotCIRaw(tday(t1)', RmE(t1)', RlE(t1)', RhE(t1)', 'r');
plot(tday(t1), ones(1, length(t1)), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$\max R_j$', '', '$R$', '', '$E$');
legend('boxoff');
ylabel('$X(t)$', 'FontSize', fnt);
xlim(tday(trange1)); ylim([0.7 2.1]); 

% Probability of resurgence
subplot(4, 2, 7);
hold on; 
plot(tdate(t1), p1Rjmax(t1), '-', 'color', 'c', 'linewidth', 2);
plot(tdate(t1), p1R(t1), '-', 'color', 'b', 'LineWidth', 2);
plot(tdate(t1), p1E(t1), '-', 'color', 'r', 'LineWidth', 2);
plot(tdate(t1), 0.5*ones(1, length(t1)), '--', 'color', 'k', 'LineWidth', 1);
h = gca; h.YColor = h.XColor;
grid off; box off; hold off;
ylabel('P($X(t) > 1$)', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt); 
xlim(tdate(trange1)); ylim([0 1.05]); h = gca;
h.XTick = [tdate(t1(1)), tdate(ceil(median(t1))), tdate(t1(end))];

% Epidemic trajectories
subplot(4, 2, 2);
plot(tdate(t2)', Ideme(:, t2)', 'g-', 'LineWidth', 2);
hold on; yymax = max(max(Ideme(:, t2)));
plot(tdate(t2)', yymax*Itot(t2)'/max(Itot(t2)), 'k-', 'LineWidth', 2);
ylabel('$I_j(t)$', 'FontSize', fnt);
grid off; box off; hold off;
xlim(tdate(trange2)); h = gca;
h.XTick = [tdate(t2(1)), tdate(ceil(median(t2))), tdate(t2(end))];

% Get ticks for middle panel
xt = h.XTick; xtlab = h.XTickLabel;
% Find time points of ticks 
tt = zeros(size(xt));
for i = 1:length(xt)-1
    tt(i) = find(tdate == xt(i));
end
tt(length(xt)) = tday(t2(end));

% Reproduction number estimates
subplot(4, 2, [4 6]);
hold on;
plotCIRaw(tday(t2)', Rmmaxj(t2)', Rlmaxj(t2)', Rhmaxj(t2)', 'c');
plotCIRaw(tday(t2)', RLam(t2)', RLaml(t2)', RLamh(t2)', 'b');
plotCIRaw(tday(t2)', RmE(t2)', RlE(t2)', RhE(t2)', 'r');
plot(tday(t2), ones(1, length(t2)), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$\max R_j$', '', '$R$', '', '$E$');
legend('boxoff');
ylabel('$X(t)$', 'FontSize', fnt);
xlim(tday(trange2)); ylim([0.7 2.1]); 

% Probability of resurgence
subplot(4, 2, 8);
hold on; 
plot(tdate(t2), p1Rjmax(t2), '-', 'color', 'c', 'linewidth', 2);
plot(tdate(t2), p1R(t2), '-', 'color', 'b', 'LineWidth', 2);
plot(tdate(t2), p1E(t2), '-', 'color', 'r', 'LineWidth', 2);
plot(tdate(t2), 0.5*ones(1, length(t2)), '--', 'color', 'k', 'LineWidth', 1);
h = gca; h.YColor = h.XColor;
grid off; box off; hold off;
ylabel('P($X(t) > 1$)', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt);
xlim(tdate(trange2)); ylim([0 1.05]); h = gca;
h.XTick = [tdate(t2(1)), tdate(ceil(median(t2))), tdate(t2(end))];
