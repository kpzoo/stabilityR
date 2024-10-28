% Bayesian predictions over horizons of chosen length
function [Rhoriz, Ihoriz, Isamp] = predHorizSamp(Rgrid, nPts, eta, horiz, qR, Lam, Iday, pR0, distvals, nday)

% Assumptions and notes
% - predicts incidence and reproduction numbers
% - returns prediction samples
% - replaces medians with means in Rhoriz and Ihoriz

%% Prediction of R across time

% Precompute state distributions
pstate = zeros(nPts, nPts);
for j = 1:nPts
    pstate(j, :) = normpdf(Rgrid(j), Rgrid, sqrt(Rgrid)*eta);
end

% Median and CIs over R
Rmed = zeros(1, horiz); Rlow = Rmed; Rhigh = Rmed;
% Distribution of R over Rgrid
Rdist = zeros(horiz, nPts); 
% Initialise with last distribution-based prediction
Rdist(1, :) = qR*pstate; Rdist(1, :) = Rdist(1, :)/sum(Rdist(1, :));

% Update this distribution over horizon
for i = 2:horiz
    Rdist(i, :) = Rdist(i-1, :)*pstate;
    % Normalise over grid
    Rdist(i, :) = Rdist(i, :)/sum(Rdist(i, :));
end

for i = 1:horiz
    % CDF and quantiles
    Rcdf = cumsum(Rdist(i, :));
    ids(1) = find(Rcdf >= 0.5, 1, 'first');
    ids(2) = find(Rcdf >= 0.025, 1, 'first');
    ids(3)= find(Rcdf >= 0.975, 1, 'first');
    clear('Rcdf'); Rmed(i) = Rgrid(ids(1)); 
    Rlow(i) = Rgrid(ids(2)); Rhigh(i) = Rgrid(ids(3));
end

% Mean and variance of distributions
Rmean = Rdist*Rgrid'; Rvar = Rdist*((Rgrid.^2)') - Rmean.^2;
% Prior variance
Rvar0 = pR0*((Rgrid.^2)') - (pR0*Rgrid').^2;
% Extract R statistics
Rstats.mean = Rmean; Rstats.var = Rvar; Rstats.var0 = Rvar0;
% Quantiles over R
Rhoriz = [Rlow' Rstats.mean Rhigh'];

%% Prediction of I across time

% Sample size for drawing incidence predictions
nsamp = 5000; Isamp = zeros(horiz, nsamp);

% First case has exact Lam and computed directly
Rsamp = datasample(Rgrid, nsamp, 'Weights', Rdist(1, :));
Isamp(1, :) = poissrnd(Rsamp*Lam(end), [1 nsamp]);

% Serial distribution over all tday
serial = serialDistrTypes(nday+horiz, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Sample total infectiousness and hence incidence
for i = 2:horiz
    % Sample from R distribution
    Rsamp = datasample(Rgrid, nsamp, 'Weights', Rdist(i, :));
    % Relevant part of serial distribution
    Pomegat = Pomega(1:(nday+i-1));
    
    % Non sampled part of total infectiousness
    Lsamp = Pomegat(i:end)*Iday(end:-1:1)'; Lsamp = Lsamp*ones(1, nsamp);
    PomegaRem = Pomegat(1:i-1);
    
    % Sample total infectiousness
    for j = 1:i-1
       Lsamp = Lsamp + PomegaRem(i-j)*Isamp(j, :);
    end
    % Project to incidence distribution
    Isamp(i, :) = poissrnd(Rsamp.*Lsamp);
end

% Statistics of predictions
Istats.mean = mean(Isamp, 2); Istats.var = var(Isamp, [], 2);
% Quantiles of incidence
Ihoriz = quantile(Isamp', [0.025, 0.975])';
Ihoriz = [Ihoriz(:, 1) Istats.mean, Ihoriz(:, end)];

