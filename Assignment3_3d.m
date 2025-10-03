data = dlmread('XYZ.txt', ',');   
X = data(:,1);    
Y = data(:,2);  
Z = data(:,3);    
N = size(data,1);

A = [X, Y, ones(N,1)];   
beta = A \ Z;             

a = beta(1);
b = beta(2);
c = beta(3);

residuals = Z - A*beta;
SSE = sum(residuals.^2);

sigma2_mle = SSE / N;        
sigma2_unbiased = SSE / (N-3); 

fprintf('Estimated plane: z = %.4f * x + %.4f * y + %.4f\n', a, b, c);
fprintf('Number of points: %d\n', N);
fprintf('MLE noise variance = %.6f\n', sigma2_mle);
fprintf('Unbiased noise variance = %.6f\n', sigma2_unbiased);

figure;
scatter3(X, Y, Z, 15, 'b', 'filled'); hold on; 
[xgrid, ygrid] = meshgrid(linspace(min(X), max(X), 30), ...
                          linspace(min(Y), max(Y), 30));
zgrid = a*xgrid + b*ygrid + c;

surf(xgrid, ygrid, zgrid, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'interp');

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Fitted Plane and Data Points');
legend('Data points', 'Fitted plane');
grid on; view(45, 30);
