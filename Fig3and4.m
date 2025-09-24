% Fig 3 and 4: stability threshold scenarios over two demes
clearvars; clc; close all; tic;

% Assumptions and notes
% - 1 deme open loop, 1 deme controlled by feedback
% - 3 scenarios with R=1 but diverse group behaviour

% Directory of some main code and plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10); 

%% Setup epidemic parameters for two demes and aggregate

% Disease of choice 
diseases = {'COVID-19', 'EBV'}; disChoice = 1;
% Setup relevant serial interval distribution and time
distvals.type = 2; nday = 201; tday = 1:nday;
switch(disChoice)
    case 0
        % COVID-19 serial interval
        distvals.pm = (1/0.65)^2; distvals.omega = 6.5;
    case 1
        % Ebola virus serial interval
        distvals.pm = 2.7066; distvals.omega = 15.3;
end
disp(['Analysing epidemics of type: ' diseases{disChoice+1}]);

% Serial distribution over all tday
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Take 95% quantiles of epidemic information
alpha = 0.025; quans = [alpha, 0.5, 1-alpha]; 

% No. runs, initial cases, R0 
M = 1; I0 = 10; R0 = 3; nDeme = 2;

%% Simulate deme 1 open loop

% Open loop Rt for deme
Rch = [1.4 1.3]+0.1; tch = 4;
R1 = Rch(1) + Rch(2)*sind(tch(1)*(1:nday));

% Daily incidence and infectiousness
I1 = zeros(1, nday); L1 = I1; 
% Initialise epidemic and warning
I1(1) = I0; Iwarn = 0; 

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    L1(i) = sum(I1(i-1:-1:1).*Pomegat);    
    % Renewal incidence
    I1(i) = poissrnd(R1(i)*L1(i));
end

% Reference signals for deme 2
ref = max(I1) - I1;

%% Control deme 2 to achieve servo target

% Outputs of controlled epidemics by deme
I2 = zeros(M, nday); L2 = I2; R2 = I2; state = I2;

% Repeat control law M times
for i = 1:M
    % Main control function
    [Is, Ls, Rs, states] = epiTracking(Pomega, nday, R0, I0, ref);
    I2(i, :) = Is;  L2(i, :) = Ls; R2(i, :) = Rs; state(i, :) = states;
end

% Sum of incidence from open loop and controlled
Isum = I1 + I2; Lsum = L1 + L2;

%% Estimate per deme and the total 

% Grid limits and noise level and grid size m
Rmin = 0.01; Rmax = 10; eta = 0.1; m = 1000; 
% Distributions (filtering and smoothing)
qR = cell(1, nDeme+1); pR = qR;
% Uniform prior over the grid of size m
p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Smoothed estimates of individual demes
[Re1, Ip1, ~, qR{1}, pR{1}] = smoothEstPred(Rgrid, m, eta, nday, p0, L1, I1);
[Re2, Ip2, ~, qR{2}, pR{2}] = smoothEstPred(Rgrid, m, eta, nday, p0, L2, I2);
% Overall estimates from sum of data
[Resum, Ipsum, ~, qR{nDeme+1}, pR{nDeme+1}] = smoothEstPred(Rgrid, m, eta,...
    nday, p0, Lsum, Isum);

% Weighted mean R based on total infectiousness
Rsum = (L1./Lsum).*R1 + (L2./Lsum).*R2;

% Posterior distribution of E and samples
qE = zeros(nday, m); pE = qE; nsamps = 20000; 

% E optimal distribution by sampling
RmE = zeros(1, nday); RlE = RmE; RhE = RmE;
for i = 1:nday
    % Individual distribution samples
    xDeme = zeros(nDeme, nsamps); xpDeme = xDeme;
    for j = 1:nDeme
        xDeme(j, :) = datasample(Rgrid, nsamps, 'Weights', qR{j}(i, :));
        xpDeme(j, :) = datasample(Rgrid, nsamps, 'Weights', pR{j}(i, :));
    end
    % E optimal samples for this day
    Dsamp = mean(xDeme); Esamp = mean(xDeme.^2)./Dsamp;
    Dpsamp = mean(xpDeme); Epsamp = mean(xpDeme.^2)./Dpsamp;

    % Statistics of E designs
    RmE(i) = mean(Esamp);
    Equants = quantile(Esamp, [0.025, 0.975]);
    RlE(i) = Equants(1); RhE(i) = Equants(2);

    % Distribution from histogram
    [qE(i, :), ~] = ksdensity(Esamp, Rgrid); qE(i, :) = qE(i, :)/sum(qE(i, :));
    [pE(i, :), ~] = ksdensity(Epsamp, Rgrid); pE(i, :) = pE(i, :)/sum(pE(i, :));
end

%% Take points of interest and compute projections

% Define horizon and times of interest
horiz = 21; tid = [50, 100, 200]; nid = length(tid);
%horiz = 21; tid = [58, 100, 196]; nid = length(tid);
% Horizon estimates and times
Rhsum = cell(1, nid); R1h = Rhsum; R2h = Rhsum; Eh = Rhsum;
Ihsum = Rhsum; I1h = Rhsum; I2h = Rhsum; IEh = Rhsum;
% Samples of projected incidence 
I1samp = Rhsum; I2samp = Rhsum; IRsamp = Rhsum; IEsamp = Rhsum;
% Times forward and backward
th = Rhsum; tb = Rhsum;

% Compute horizon estimates and predictions
for i = 1:nid
    % Times for horizon and backwards
    th{i} = tday(tid(i))+1:tday(tid(i))+horiz;
    tb{i} = tday(tid(i))-2*horiz:tday(tid(i));

    % Projections and confidence intervals
    [R1h{i}, I1h{i}, I1samp{i}] = predHorizSamp(Rgrid, m, eta, horiz,...
        pR{1}(tid(i), :), L1, I1, p0, distvals, nday);
    [R2h{i}, I2h{i}, I2samp{i}] = predHorizSamp(Rgrid, m, eta, horiz,...
        pR{2}(tid(i), :), L2, I2, p0, distvals, nday);
    [Rhsum{i}, Ihsum{i}, IRsamp{i}] = predHorizSamp(Rgrid, m, eta,...
        horiz, pR{nDeme+1}(tid(i), :), Lsum, Isum, p0, distvals, nday);
    [Eh{i}, IEh{i}, IEsamp{i}] = predHorizSamp(Rgrid, m, eta, horiz,...
        pE(tid(i), :), Lsum, Isum, p0, distvals, nday);
end

% Summing incidence samples from demes
Itotsamp = cell(1, nid); Idemeh = Itotsamp;
for i = 1:nid
    % Construct projections from sum of deme projected incidence
    Itotsamp{i} = I1samp{i} + I2samp{i};
    % Quantiles and mean
    Idemequant = quantile(Itotsamp{i}', [0.025, 0.975])';
    Idememean = mean(Itotsamp{i}, 2);
    Idemeh{i} = [Idemequant(:, 1) Idememean, Idemequant(:, end)];
end

% Horizons to consider
horid = [7 14 21]; nh = length(horid);
% Last weekly mean of infections (observed)
Ilast = mean(Isum(end:-1:end-3)); nsampI = length(Itotsamp{1}');

% Probabilities that infections are growing
IRprob = zeros(nid, horiz); IEprob = IRprob; Itotprob = IRprob;
for i = 1:nid
    Itottemp = cell(1, horiz); IRtemp = Itottemp; IEtemp = Itottemp;
    for j = 1:horiz
        % Probability of growth in incidence
        Itotprob(i, j) = sum(Itotsamp{i}(j, :) > Ilast)/nsampI;
        IRprob(i, j) = sum(IRsamp{i}(j, :) > Ilast)/nsampI;
        IEprob(i, j) = sum(IEsamp{i}(j, :) > Ilast)/nsampI;
    end
end

% Probabilities that infections are growing with samples
IRden = cell(1, nid); IEden = IRden; Itotden = IRden;
IRproj = IRden; IEproj = IRproj; Itotproj = IRproj;
for i = 1:nid
    [IRden{i}, IEden{i}, Itotden{i}] = getProjGrowthProb(Ilast, horiz, tid(i), nsampI,...
        p0, L1, I1, I2, L2, distvals, nday, eta, m, Rgrid, Lsum, Isum, pR, pE);
    % Statistics of growth projections
    Itotq = quantile(Itotden{i}', [0.025 0.975])'; Itotm = mean(Itotden{i}, 2);
    Itotproj{i} = [Itotq(:, 1) Itotm Itotq(:, 2)];
    IRq = quantile(IRden{i}', [0.025 0.975])'; IRm = mean(IRden{i}, 2);
    IRproj{i} = [IRq(:, 1) IRm IRq(:, 2)];
    IEq = quantile(IEden{i}', [0.025 0.975])'; IEm = mean(IEden{i}, 2);
    IEproj{i} = [IEq(:, 1) IEm IEq(:, 2)];
end

%% Publishable figure 3

figure('Position', [10 10 800 800]);
subplot(2, 1, 1); hold on;
plotCIRaw(tday(2:end)', Ip1.mean, Ip1.low, Ip1.high, 'g');
plotCIRaw(tday(2:end)', Ip2.mean, Ip2.low, Ip2.high, 'g');
plotCIRaw(tday(2:end)', Ipsum.mean, Ipsum.low, Ipsum.high, 'b');
xlim([tday(2) tday(end)]); h = gca; ylimh = h.YLim;
for i = 1:length(tid)
    plot([tid(i) tid(i)], ylimh, 'k--', 'LineWidth', 2);
end
plot(tday, I1, 'gx', 'MarkerSize', 5, 'LineWidth', 1);
plot(tday, I2, 'gx', 'MarkerSize', 5, 'LineWidth', 1);
plot(tday, I1+I2, 'bx', 'MarkerSize', 5, 'LineWidth', 1);
hold off; box off; grid off; xlim([21 201]);
ylabel('predictions $I(t)$', 'FontSize', fnt);
title('One-step-ahead incidence', FontSize = fnt);

subplot(2, 1, 2); hold on;
plotCIRaw(tday', Re1.mean, Re1.low, Re1.high, 'g');
plotCIRaw(tday', Re2.mean, Re2.low, Re2.high, 'g');
plotCIRaw(tday', Resum.mean, Resum.low, Resum.high, 'b');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
xlim([tday(2) tday(end)]); h = gca; ylimh = h.YLim;
for i = 1:length(tid)
    plot([tid(i) tid(i)], ylimh, 'k--', 'LineWidth', 2);
end
plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2);
hold off; box off; grid off; 
xlim([21, 201]); ylim([0, 4]);
ylabel('estimates $X(t)$', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt);
title('Reproduction numbers', FontSize = fnt);
legend('$R_j$', '', '', '', '$R$', '', 'E', 'Location',...
    'best', fontsize = fnt);

%% Publishable figure 4

% Distributions and predictions
figure('Position', [10 10 1000 800]);
ylabs = {'$\textbf{A} \,$', '$\textbf{B} \,$', '$\textbf{C} \,$'};
mark = 5; smth = 5;
for i = 1:nid
    
    % Past incidence (remove +1 as 1-step-ahead at that time)
    subplot(nid, 3, 3*(i-1)+1);
    plot(tb{i}, smooth(Isum(tb{i}), smth), 'b', 'LineWidth', 2);
    hold on;
    plot(tb{i}, Isum(tb{i}), 'bx', 'MarkerSize', mark, 'LineWidth', 1);
    plot(tb{i}, smooth(I1(tb{i}), smth), 'g', 'LineWidth', 2);
    plot(tb{i}, smooth(I2(tb{i}), smth), 'g', 'LineWidth', 2);
    plot(tb{i}, I1(tb{i}), 'gx', 'MarkerSize', mark, 'LineWidth', 1);
     plot(tb{i}, I2(tb{i}), 'gx', 'MarkerSize', mark, 'LineWidth', 1);
    h = gca; plot(tid(i)*ones(1, 2), h.YLim, 'k--', 'LineWidth', 2);
    hold off; box off; grid off;
    xlim([tb{i}(1) tb{i}(end)]);
    if i == 1
        title('incidence $I_1^{\tau}$', 'FontSize', fnt);
    end
    if i == nid
        xlabel('$t$ (days)', 'FontSize', fnt);
    end
    ylabel(ylabs{i}, 'FontSize', fnt, 'Rotation', 0);


    % Histograms of filtering distributions
    subplot(nid, 3, 3*(i-1)+2);
    plot(Rgrid, pR{1}(tid(i), :), 'Color', 'g', 'LineWidth', 2);
    hold on;
    plot(Rgrid, pR{2}(tid(i), :), 'Color', 'g', 'LineWidth', 2);
    plot(Rgrid, pR{nDeme+1}(tid(i), :), 'Color', 'b', 'LineWidth', 2);
    plot(Rgrid, pE(tid(i), :), 'Color', 'r', 'LineWidth', 2);
    h = gca;
    plot([1 1], h.YLim, 'k--', 'LineWidth', 1);
    hold off; box off; grid off;
    xlim([0 4]);
    if i == 1
        title('posterior P$(X | I_1^{\tau})$', 'FontSize', fnt);
    end
    if i == nid
        xlabel('$X$ at time $\tau$', 'FontSize', fnt);
        legend('$R_1$', '$R_2$', '$R$', '$E$', 'Location','Best');
        legend boxoff;
    end

    % Probability of growth with time
    subplot(nid, 3, 3*(i-1)+3);
    plotCIRaw(th{i}', Itotproj{i}(:, 2), Itotproj{i}(:, 1), Itotproj{i}(:, 3), 'g');
    hold on;
    plotCIRaw(th{i}', IRproj{i}(:, 2), IRproj{i}(:, 1), IRproj{i}(:, 3), 'b');
    plotCIRaw(th{i}', IEproj{i}(:, 2), IEproj{i}(:, 1), IEproj{i}(:, 3), 'r');
    
    hold off; box off; grid off;
    xlim([th{i}(1) th{i}(end)]);
    if i == 1
        title('P($I(\tau+h) > I(\tau)$)', 'FontSize', fnt);
    end
    if i == nid
        xlabel('$\tau+h$ (days)', 'FontSize', fnt);
    end
end


% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
