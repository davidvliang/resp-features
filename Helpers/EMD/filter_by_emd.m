function [signal, K_optim, imf, info] = filter_by_emd(Bpha_uw)
%FILTER_BY_EMD Summary of this function goes here
%   Detailed explanation goes here

% Compute EMD and obtain IMF
[imf, residual, info] = emd(Bpha_uw);
% emd(Bpha_uw);

% Calculate Mutual Information MI(k) via Fast MI
MI = zeros(size(imf,2)-1,1);
Xr = zeros(size(imf,1),1); % deterministic component (respiratory)
Xn = zeros(size(imf,1),1); % stochastic component (noise)

for idx = 1:size(MI,1)
    K_temp = idx + 1;
    Xr = sum(imf(:, (K_temp:size(imf,2))), 2)+residual; % add imfs k through m
    Xn = sum(imf(:, (1:K_temp-1)), 2)+residual; % add imfs 1 through k-1
%     MI(idx) = mi(Xr,Xn); % Fast MI
    MI(idx) = mi_cont_cont(Xn, Xr, 1); % knn method
end

% Calculate Mutual Information Ratio MIR(k) [eq. 9]
MIR = zeros(size(MI,1)-1,1);
for idx = 1:size(MI,1)-1
    MIR(idx) = MI(idx+1) / MI(idx);
end

% Find optimal K value (w/ highest MIR)
[~, mir_argmax] = max(MIR);
K_optim = mir_argmax+1;

% Reconstruct the filtered signal [eq. 6]
signal = sum(imf(:, (K_optim:size(imf,2))), 2)+residual;

end

