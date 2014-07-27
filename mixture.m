% Created by: Sherif Abdelwahab
% Updated: Nov 20, 2013
% Simulate the mixture model

function [kpq]=mixture(n,N,p,q,a)

kp = binornd(n,p,N,1); % generate N samples from binomial dist with p
kq = binornd(n,q,N,1); % generate N samples from binomial dist with q
ind = [rand(N,1)<a];  % indicator to select either kp or kq with prob a
kpq = ind .* kp + (1 - ind) .* kq;