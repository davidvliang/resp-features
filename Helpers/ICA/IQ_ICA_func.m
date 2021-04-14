% Uses JADE to find source signals from two sources in IQ modulated form
% relies on pca_lin.m to demodulate IQ signals
% and jader.m to perform JADE ICA

% in: column vectors of equal length
% out: normalized 2xT matrix representing 2 independent sources
function [source1, source2] = IQ_ICA_func(I1, Q1, I2, Q2)

% Demodulate the IQ signals using Mahmud's code
IQ1_demod = pca_lin(I1, Q1);
IQ2_demod = pca_lin(I2, Q2);

% format the imput correctly for the jade function
X = [IQ1_demod'; IQ2_demod'];
B = jader(X);

% compute the results
S = B*X;
source1 = -normalize(S(1,:)');
source2 = -normalize(S(2,:)');
