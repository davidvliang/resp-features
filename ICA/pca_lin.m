function x_lin = pca_lin (I, Q);
% linear demodulation
% I and Q should be column vectors
x = [I Q];
[n p] = size (x);
x = x - ones (n,1)* mean(x); % DC removal
cov_x = x'*x/(n-1);
[v D] = eig (cov_x);
x = x * v (:,2);
x_lin = x;