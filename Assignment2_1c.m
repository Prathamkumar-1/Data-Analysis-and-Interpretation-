clear; close all; clc;

% fit_plane_XYZ.m
% Reads XYZ.txt using dlmread (comma-separated), 
% fits z = a*x + b*y + c, and prints the results.

% --- Read data ---
data = dlmread('XYZ.txt', ',');   % specify comma delimiter
X = data(:,1);    % first column
Y = data(:,2);    % second column
Z = data(:,3);    % third column
N = size(data,1);

% --- Build design matrix and solve ---
A = [X, Y, ones(N,1)];    % N x 3
beta = A \ Z;             % least squares solution [a; b; c]

a = beta(1);
b = beta(2);
c = beta(3);

% --- Residuals and noise variance ---
residuals = Z - A*beta;
SSE = sum(residuals.^2);

sigma2_mle = SSE / N;          % MLE variance
sigma2_unbiased = SSE / (N-3); % unbiased estimate

% --- Print results ---
fprintf('Estimated plane: z = %.4f * x + %.4f * y + %.4f\n', a, b, c);
fprintf('Number of points: %d\n', N);
fprintf('MLE noise variance = %.6f\n', sigma2_mle);
fprintf('Unbiased noise variance = %.6f\n', sigma2_unbiased);
