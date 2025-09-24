% Fig 2: square error convergence of global statistics
clearvars; clc; close all; 

% Assumptions and notes
% - sample EpiEstim gamma posteriors to get Rj
% - sample different Dirichlet weightings
% - compute MSE for Rj combined, Rmax, R and E

% Directory and where saving
thisDir = cd; saveFol = 'results/';
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);


%% Simulate estimated reproduction numbers

% Samples and number of demes 
ps = [5 10 20 40 60 80]; 
np = length(ps); nsamps = 5000; 

% Weights of Lj variations and num runs
wvar = 4; nw = 1000; 
% Estimates of R with mean 1
Ltot = 1000; Itot = Ltot;

% Performance indices
Rsq = zeros(nw, np); Esq = Rsq; 
Rjsq = Rsq; Rmaxsq = Rsq; idorder = Rsq; 
cE = zeros(1, np); cmax = cE; 

% Performance with p
for k = 1:np
    p = ps(k);

    % Equally divided Lj and Dirichlet weights for Ij
    Lj = Ltot/p; pw = wvar*ones(1, p);
    % Random draws of weights summing to 1
    wj = drchrnd(pw, nw); Rjmeantheo = wj*p;

    % For every weight set sample Rj and also find E
    Rjset = cell(1, nw); E = zeros(nw, nsamps);
    Rjmean = zeros(nw, p); Rsamp = E; Rmax= E;
    for i = 1:nw
        % New set of Rj possibilities
        Rj = zeros(p, nsamps);
        for j = 1:p
            Rj(j, :) = gamrnd(wj(i, j)*Itot, 1/Lj, [1 nsamps]);
        end
        % E number from the Rj
        E(i, :) = sum(Rj.*Rj)./sum(Rj);
        % Max Rj statistic
        Rmax(i, :) = max(Rj);
        % Save individual Rj samples and sampled mean
        Rjset{i} = Rj; Rjmean(i, :) = mean(Rjset{i}, 2);
        % Recompute R distribution from Rj
        Rsamp(i, :) = sum(Lj*Rj/Ltot);
    end

    % Order by how different Rj are from R
    Rjsq(:, k) = mean((Rjmean - 1).^2, 2); 
    [~, idorder(:, k)] = sort(Rjsq(:, k));

    % Square error for all sample sets (run at nw = 100)
    Rsq(:, k) = mean((Rsamp - 1).^2, 2); 
    Esq(:, k) = mean((E - 1).^2, 2);
    Rmaxsq(:, k) = mean((Rmax - 1).^2, 2);
    % Correlation of E or Rmax and Rj mean square errors
    cE(k) = corr(Esq(:, k), Rjsq(:, k));
    cmax(k) = corr(Rmaxsq(:, k), Rjsq(:, k));
    cE(k) = roundn(cE(k), -2); cmax(k) = roundn(cmax(k), -2);
end

%% Publishable figure

% Metrics across order and runs
figure('Position', [10 10 800 800]);
for k = 1:np
    subplot(np/2, 2, k);
    plot(1:nw, Rjsq(idorder(:, k), k), 'g', 'LineWidth', 2);
    hold on;
    plot(1:nw, Rsq(idorder(:, k), k), 'b', 'LineWidth', 2);
    plot(1:nw, Esq(idorder(:, k), k), 'r', 'LineWidth', 2);
    %plot(1:nw, Rmaxsq(idorder(:, k), k), 'c', 'LineWidth', 2);
    hold off; grid off; box off;
    xlabel(['runs $| \, p = $ ' num2str(ps(k))], 'FontSize', fnt);
    ylabel(['$\rho = $ ' num2str(cE(k)) ', ' num2str(cmax(k))], 'FontSize', fnt);
    if k == 1 || k == 2
        title('MSE from 1', 'FontSize', fnt);
        legend('$R_j$', '$R$', '$E$', 'Location', 'best', fontsize = fnt);
    end
end









