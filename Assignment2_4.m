clear; close all; clc;

%% Reading the images and storing them as 2D arrays (and casting into double)
I1 = double(imread('C:\Users\Pratham Kumar\Downloads\T1.jpg'));
I2_orig = double(imread('C:\Users\Pratham Kumar\Downloads\T2.jpg'));
[H, W] = size(I1); %H is height and W is weight

%% Fixed parameters
txRange = -10:10;
binWidth = 10;
nbins = ceil(256/binWidth);

%% All three cases of I2
I2_A = I2_orig;
I2_B = 255 - I1;
I2_C = 255*(I1.^2)/max(I1(:).^2);

%% Initializing result arrays
rhoA = nan(size(txRange)); %correlation coeff
qmiA = nan(size(txRange)); %QMI 
miA = nan(size(txRange));  %MI

rhoB = nan(size(txRange)); 
qmiB = nan(size(txRange));
miB = nan(size(txRange));

rhoC = nan(size(txRange));
qmiC = nan(size(txRange));
miC = nan(size(txRange));

for idx = 1:3
    switch idx
        case 1
            I2base = I2_A;
            rhoVec = rhoA;
            qmiVec = qmiA;
            miVec = miA;
        case 2
            I2base = I2_B;
            rhoVec = rhoB;
            qmiVec = qmiB;
            miVec = miB;
        case 3
            I2base = I2_C;
            rhoVec = rhoC;
            qmiVec = qmiC;
            miVec = miC;
    end
    
    %shifting I2 image by tx units
    for a = 1:numel(txRange)
        tx = txRange(a);
        I2s = zeros(H, W);
        mask = false(H, W); %intialize mask to be all false
        if tx >= 0
            I2s(:, 1+tx:W) = I2base(:, 1:W-tx);
            mask(:, 1+tx:W) = true;  %set these shifted points as valid(mask = true).
        else
            sh = -tx;
            I2s(:, 1:W-sh) = I2base(:, 1+sh:W);
            mask(:, 1:W-sh) = true; %mark the shifted pixels as valid(mask = true).
        end
        x = I1(mask);
        y = I2s(mask);
        
        % Computing correlation coeff
        if isempty(x)
            rhoVal = NaN;
        else
            xm = mean(x);
            ym = mean(y);
            num = sum((x - xm).*(y - ym));
            den = sqrt(sum((x - xm).^2) * sum((y - ym).^2));
            if den == 0    
                  rhoVal = NaN;
            else
                rhoVal = num/den;
            end

        end
        
        P12 = zeros(nbins, nbins);
        xb = floor(x./binWidth) + 1;
        yb = floor(y./binWidth) + 1;
        xb(xb<1) = 1; 
        xb(xb>nbins) = nbins;
        yb(yb<1) = 1; 
        yb(yb>nbins) = nbins;
        
        for k = 1:numel(xb)
            P12(xb(k), yb(k)) = P12(xb(k), yb(k)) + 1;
        end
        
        s = sum(P12(:));
        if s > 0
            P12 = P12 / s;
        end
        
        P1 = sum(P12, 2);
        P2 = sum(P12, 1);
        P1P2 = P1 * P2;
        
        % Computing QMI
        QMIval = sum((P12(:) - P1P2(:)).^2);
        
        % Computing MI
        nz = P12 > 0;
        MIval = sum(P12(nz) .* log2(P12(nz) ./ P1P2(nz)));
        
        switch idx
            case 1
                rhoA(a) = rhoVal;
                qmiA(a) = QMIval;
                miA(a) = MIval;
            case 2
                rhoB(a) = rhoVal;
                qmiB(a) = QMIval;
                miB(a) = MIval;
            case 3
                rhoC(a) = rhoVal;
                qmiC(a) = QMIval;
                miC(a) = MIval;
        end
    end
end

%% Plotting
figure; 
plot(txRange, rhoA, '-o', txRange, rhoB, '-s', txRange, rhoC, '-d');
grid on;
xlabel('tx (pixels)');
ylabel('\rho');
title('Correlation vs tx');
legend('Case A: I2 = T2.jpg', 'Case B: I2 = 255 - I1', 'Case C: nonlinear', 'Location', 'best');

figure;
plot(txRange, qmiA, '-o', txRange, qmiB, '-s', txRange, qmiC, '-d');
grid on;
xlabel('tx (pixels)');
ylabel('QMI');
title('QMI vs tx');
legend('Case A', 'Case B', 'Case C', 'Location', 'best');

figure;
plot(txRange, miA, '-o', txRange, miB, '-s', txRange, miC, '-d');
grid on;
xlabel('tx (pixels)');
ylabel('MI (bits)');
title('Mutual Information vs tx');
legend('Case A', 'Case B', 'Case C', 'Location', 'best');

