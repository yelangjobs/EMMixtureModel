% Created by: Sherif Abdelwahab
% Updated: Nov 23, 2013
% Computing the MOM estimator using the factorial moment method

function [theta_hat,error,exec_time]=mom(theta,n,N,theta_g,kpq)

tic;

a0 = theta(1); % exact value alpha
p0 = theta(2); % exact value p
q0 = theta(3); % exact value q



%% MOM algorithm 

% generate N data sequence of binomial mixture model

% kp = binornd(n,p0,N,1)'; % generate N samples from binomial dist with p
% kq = binornd(n,q0,N,1)'; % generate N samples from binomial dist with q
% 
% ind = [rand(N,1)'<a0];  % indicator to select either kp or kq with prob a
% kpq = ind .* kp + (1 - ind) .* kq; %mixture model at value alpha=a



[row,col]=size(kpq);

M1 = [];
M2 = [];
M3 = [];
for i=[1:row]
    k = kpq(i,1);
     
    M1 = [M1 k ];
    M2 = [M2 k * (k-1) ];
    M3 = [M3 k * (k-1) * (k-2) ];
       
    
end

M1 = sum(M1) / N / n;
M2 = sum(M2) / N / n / (n-1);
M3 = sum(M3) / N / n / (n-1) / (n-3);

Z = (M3 - M2 * M1) / (M2 - M1^2);
quantity = Z^2-4*M1*Z+4*M2;
root = sqrt(Z^2-4*M1*Z+4*M2);


if quantity > 0 && root <= min(Z,2-Z)
    if q0>=p0
    pe = 0.5 * (Z - root);

    qe = 0.5 * (Z + root);

    ae = (M1 - qe)/(pe - qe);
    else
    pe = 0.5 * (Z + root);

    qe = 0.5 * (Z - root);

    ae = 1-(M1 - qe)/(pe - qe);
    end
else
    pe = M1;
    qe = M1;
    ae = 0;
end

theta_hat=[ae;pe;qe];
error=theta_hat-theta;

exec_time=toc;


