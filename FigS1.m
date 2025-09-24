% Fig S1: multi-deme E and R with empirical data from Texas, USA
clearvars; clc; close all; tic;

% Assumptions and notes
% - accounts for different starting times
% - optimal designs achieved by sampling from qR distributions
% - no interactiion among demes and same serial intervals used
% - data from https://github.com/nytimes/covid-19-data

% Directory and where saving (or loading)
thisDir = cd; saveFol = 'data/texas'; 
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Load incidence data and serial intervals
cd(saveFol);
% Read state data
files = dir("Texas2020.csv");
data = readtable(files.name);
cd(thisDir);

% Get cases from demes comprising 80% of all infections
[tdate, Ideme, nDeme, nFull] = procUSdata(data, 0.75);
nday = length(tdate); tday = 1:nday;

%% Format empirical data and obtain total infectiousness

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

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

%% Figures on overall and individual deme performance

% Time to narrow down simulations on
tstart = 1; tstop = nday;
trange = [tstart tstop];

% Peak sizes of epidemics
figure('Position', [10 10 800 1000]);
subplot(3, 1, 1);
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
subplot(3, 1, [2 3]);
hold on;
%plotCIRaw(tday', Rmmaxj', Rlmaxj', Rhmaxj', 'c');
plotCIRaw(tday', RLam', RLaml', RLamh', 'b');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$R$', '', '$E$');
legend('boxoff');
ylabel('$X(t)$', 'FontSize', fnt);
xlim(tday(trange)); 
xlabel('$t$ (days)', 'FontSize', fnt);

%% Publishable Figure S1

% Cut points for investigation
tstart = find(tdate == '2020-04-01');
tstop = find(tdate == '2020-06-10');
% Corresponding ranges for investigation
trange1 = [tstart tstop]; t1 = tstart:tstop;

figure('Position', [10 10 1000 800]);
% Epidemic trajectories
subplot(4, 1, 1);
% Total cases per deme in the window
demeTotals = sum(Ideme(:, t1),2);
% Sort demes by total descending
[~, sortIdx] = sort(demeTotals,'descend');
IdemeSorted  = Ideme(sortIdx, t1);
area(tdate(t1), IdemeSorted'); hold on;
plot(tdate(t1), Itot(t1), 'k-', 'LineWidth', 2);
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
subplot(4, 1, [2 3]);
hold on;
plotCIRaw(tday(t1)', RLam(t1)', RLaml(t1)', RLamh(t1)', 'b');
plotCIRaw(tday(t1)', RmE(t1)', RlE(t1)', RhE(t1)', 'r');
plot(tday(t1), ones(1, length(t1)), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$R$', '', '$E$');
legend('boxoff');
ylabel('$X(t)$', 'FontSize', fnt);
xlim(tday(trange1)); ylim([0.7 2.1]); 

% Probability of resurgence
subplot(4, 1, 4);
hold on; 
plot(tdate(t1), p1R(t1), '-', 'color', 'b', 'LineWidth', 2);
plot(tdate(t1), p1E(t1), '-', 'color', 'r', 'LineWidth', 2);
plot(tdate(t1), 0.5*ones(1, length(t1)), '--', 'color', 'k', 'LineWidth', 1);
h = gca; h.YColor = h.XColor;
grid off; box off; hold off;
ylabel('P($X(t) > 1$)', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt); 
xlim(tdate(trange1)); ylim([0 1.05]); h = gca;
h.XTick = [tdate(t1(1)), tdate(ceil(median(t1))), tdate(t1(end))];
