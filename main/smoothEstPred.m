function [Rest, Ipred, prL1S, qR, pR] = smoothEstPred(Rgrid, m, eta, nday, p0, Lam, Iloc)

% Assumptions and notes
% - added output of filtering pR
% - runs filter and smoother in one go
% - includes growth rates with confidence

% EpiFilter estimates using local cases
[~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Lam, Iloc);
% EpiSmoother estimates for single trajectory
[~, RlowS, RhighS, RmeanS, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiSmoother one-step-ahead predictions 
[predS, predIntS] = recursPredict(Rgrid, qR, Lam, RmeanS, max(Iloc));
predlowS = predIntS(:,1)'; predhighS = predIntS(:,2)';

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last'); prL1S = zeros(1, nday); 
% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and prob R <= 1
    Rcdf = cumsum(qR(i, :)); prL1S(i) = Rcdf(id1);
end

% Output data structures for R estimates
Rest.mean = RmeanS'; Rest.low = RlowS'; Rest.high = RhighS';
% Output data structures for I predictions
Ipred.mean = predS'; Ipred.low = predlowS'; Ipred.high = predhighS';