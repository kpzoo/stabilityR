% Fig 1A analysis 
clearvars; clc; close all;

% Dimension of range
nR = 1000; nw = 10; 
% Global statistics
R = zeros(nw, nR); E = R; Rmax = R;

% Local statistics
R1 = linspace(0, 2, nR); R2 = 2 - R1;
w1 = linspace(0, 1/2, nw); w2 = 1 - w1; 

% Compute R and E from local values
for i = 1:nw
    for j = 1:nR
        R(i, j) = w1(i)*R1(j) + w2(i)*R2(j);
        E(i, j) = (R1(j)^2 + R2(j)^2)/(R1(j) + R2(j));
        Rmax(i, j) = max(R1(j), R2(j));
    end
end

figure;
plot(R1, R(end, :), 'b', 'LineWidth', 2);
hold on;
plot(R1, E(end, :), 'r', 'LineWidth', 2);
plot(R1, Rmax(end, :), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
ylim([0.9 2]); xlim([0.8 2]);
xlabel('$R_1$'); ylabel('$X$');

