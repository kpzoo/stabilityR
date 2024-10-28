% Fig 1C analysis: gamma distributions for R = 1 and E = 1
clearvars; clc; close all; tic;

% Assumptions and notes
% - optimal designs achieved by sampling EpiEstim gamma posteriors
% - examine at two scenarios of L1 = 20 and 100

% Directory and where saving
thisDir = cd; saveFol = 'Results/';
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);


%% Simulate estimated reproduction numbers

% Samples and number of demes 
nsamps = 5000; p = 5; 

% Weights of Lj variations and num runs
varyRj = 1; nw = 4; 

% Estimates of R with mean 1
Ltot = 200; Itot = Ltot;
R = gamrnd(Itot, 1/Ltot, [1 nsamps]);

% Different variations on weigths of Lj
switch(varyRj)
    case 0
        % Low variation
        wvar = 200;
    case 1
        % High variation
        wvar = 4;
end

% Equally divided Lj and Dirichlet weights for Ij
Lj = Ltot/p; pw = wvar*ones(1, p); 
% Random draws of weights summing to 1
wj = drchrnd(pw, nw); Rjmeantheo = wj*p;

% For every weight set sample Rj and also find E
Rjset = cell(1, nw); E = zeros(nw, nsamps); 
Rjmean = zeros(nw, p); Rsamp = E;
for i = 1:nw
    % New set of Rj possibilities
    Rj = zeros(p, nsamps);
    for j = 1:p
        Rj(j, :) = gamrnd(wj(i, j)*Itot, 1/Lj, [1 nsamps]);
    end
    % E number from the Rj
    E(i, :) = sum(Rj.*Rj)./sum(Rj);
    % Save individual Rj samples and sampled mean
    Rjset{i} = Rj; Rjmean(i, :) = mean(Rjset{i}, 2);
    % Recompute R distribution from Rj
    Rsamp(i, :) = sum(Lj*Rj/Ltot);
end

% Order by how different Rj are from R
Rsumsq = mean((Rjmean - 1).^2, 2); [~, idorder] = sort(Rsumsq);

if nw <= 20
    % Plot histograms of estimates
    figure('Position', [10 10 800 800]); ax = zeros(1, nw);
    for i = 1:nw
        ax(i) = subplot(nw/2, 2, i);
        hold on;
        h = histogram(R, 'BinWidth', 0.05, 'Normalization','probability');
        h.EdgeAlpha = 0; h.FaceColor = 'b'; h.FaceAlpha = 0.3;
        h = histogram(E(idorder(i), :), 'BinWidth', 0.05, 'Normalization','probability');
        h.EdgeAlpha = 0; h.FaceColor = 'r'; h.FaceAlpha = 0.3;
        h = gca; yval = h.YLim;
        plot(repmat(Rjmean(idorder(i), :), [2 1]), h.YLim, 'g--', 'LineWidth', 2);
        hold off; grid off; box off;
        if ismember(i, [nw nw-1])
            xlabel('$X$ at time $t$', 'FontSize', fnt);
        end
        if ismember(i, 2*(1:nw/2) - 1)
            ylabel('P$(X \, | \, I_1^t)$', 'FontSize', fnt);
        end
    end
    linkaxes(ax, 'xy');
end

% Version with complete histograms for all groups
if nw <= 20
    figure('Position', [10 10 800 800]); ax = zeros(1, nw);
    for i = 1:nw
        ax(i) = subplot(nw/2, 2, i);
        hold on;
        for j = 1:p
            h = histogram(Rjset{idorder(i)}(j, :), 'BinWidth', 0.05,...
                'Normalization','probability');
            h.EdgeAlpha = 0; h.FaceColor = 'g'; h.FaceAlpha = 0.1;
        end
        h = histogram(R, 'BinWidth', 0.05, 'Normalization','probability');
        h.EdgeAlpha = 0; h.FaceColor = 'b'; h.FaceAlpha = 0.3;
        h = histogram(E(idorder(i), :), 'BinWidth', 0.05, 'Normalization','probability');
        h.EdgeAlpha = 0; h.FaceColor = 'r'; h.FaceAlpha = 0.3;
        h = gca; yval = h.YLim;
        plot(repmat(Rjmean(idorder(i), :), [2 1]), h.YLim, 'g--', 'LineWidth', 1);
        hold off; grid off; box off;
        if ismember(i, [nw nw-1])
            xlabel('$X$ at time $t$', 'FontSize', fnt);
        end
        if ismember(i, 2*(1:nw/2) - 1)
            ylabel('P$(X \, | \, I_1^t)$', 'FontSize', fnt);
        end
    end
    linkaxes(ax, 'xy');
end






