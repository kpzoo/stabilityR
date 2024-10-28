% Function to sample from gamma posteriors and construct E and R
function [ms, vs, Fs, Rdeme] = getERgamma(aj, L1, ndel, lena)

% Delta to incidence on region 1 and sample size
del = linspace(0.1, 0.9, ndel); nSamp = 20000;

% Variables for comparison
m1 = zeros(lena, ndel); mR = m1; mE = m1;
v1 = zeros(lena, ndel); vR = v1; vE = v1;
F1 = zeros(lena, ndel); FR = F1; FE = F1;
Rdeme = zeros(lena, ndel);

for j = 1:lena
    % Total infectiousness in each region scaled by a
    a = aj(j); L2 = a*L1; I2 = L2;
    % Weights for the regions based on a
    w1 = L1/(L1 + L2); w2 = 1 - w1; I1 = (1 + del)*L1; 

    % Samples of R1, R (combined) and E
    R1 = zeros(ndel, nSamp); R = R1; E = R1;

    % Fix sample from region 2 as no perturbations there
    R2 = gamrnd(I2, 1/L2, [1 nSamp]);
    % Ratio of means (relative Rj number)
    Rdeme(j, :) = (I1/L1)/(I2/L2);

    % Obtain samples and estimates of consensus
    for i = 1:ndel
        % Region 1 reproduction number
        R1(i, :) = gamrnd(I1(i), 1/L1, [1 nSamp]);

        % Overall R and E number
        R(i, :) = w1*R1(i, :) + w2*R2;
        E(i, :) = (R1(i, :).^2 + R2.^2)./(R1(i, :) + R2);

        % Prob > 1 (resurgence) for all metrics
        F1(j, i) = length(find(R1(i, :) > 1))/nSamp;
        FR(j, i) = length(find(R(i, :) > 1))/nSamp;
        FE(j, i) = length(find(E(i, :) > 1))/nSamp;
    end

    % Statistics from samples
    m1(j, :) = mean(R1, 2); v1(j, :) = var(R1, [] , 2)';
    mR(j, :) = mean(R, 2); vR(j, :) = var(R, [] , 2)';
    mE(j, :) = mean(E, 2); vE(j, :) = var(E, [] , 2)';
end

% Outputs of statistics
ms.R1 = m1; ms.R = mR; ms.E = mE;
vs.R1 = v1; vs.R = vR; vs.E = vE;
Fs.R1 = F1; Fs.R = FR; Fs.E = FE;