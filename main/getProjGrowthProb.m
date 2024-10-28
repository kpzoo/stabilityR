% Determine probability of growth from projections
function [IR, IE, Itot] = getProjGrowthProb(Ilast, horiz, tact, nsampI,...
    p0, L1, I1, I2, L2, distvals, nday, eta, m, Rgrid, Lsum, Isum, pR, pE)

% Assumptions and notes
% - compare samples of projections to Ilast for E, R and sum 
% - provides mean and 95% credible intervals
% - only works with 2 demes

% Replicate sample draws
nRep = 100; nDeme = 2;
Itot = zeros(horiz, nRep); IR = Itot; IE = Itot;

% Compute horizon estimates and predictions
for i = 1:nRep
    % Local deme incidence projections (sum)
    [~, ~, I1samp] = predHorizSamp(Rgrid, m, eta, horiz,...
        pR{1}(tact, :), L1, I1, p0, distvals, nday);
    [~, ~, I2samp] = predHorizSamp(Rgrid, m, eta, horiz,...
        pR{2}(tact, :), L2, I2, p0, distvals, nday);
    Itotsamp = I2samp + I1samp;

    % Global statistic based incidence projections
    [~, ~, IRsamp] = predHorizSamp(Rgrid, m, eta,...
        horiz, pR{nDeme+1}(tact, :), Lsum, Isum, p0, distvals, nday);
    [~, ~, IEsamp] = predHorizSamp(Rgrid, m, eta, horiz,...
        pE(tact, :), Lsum, Isum, p0, distvals, nday);

    for j = 1:horiz
        % Probability of growth in incidence
        Itot(j, i) = sum(Itotsamp(j, :) > Ilast)/nsampI;
        IR(j, i) = sum(IRsamp(j, :) > Ilast)/nsampI;
        IE(j, i) = sum(IEsamp(j, :) > Ilast)/nsampI;
    end
end

