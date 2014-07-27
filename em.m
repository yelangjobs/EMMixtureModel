% Created by: Sherif Abdelwahab
% Updated: Nov 29, 2013
% Updated: Dec 11, 2013
% Implemntation of the EM algorithm for mixture of two binomial
% q>p
function [theta_hat,error,exec_time]=em(theta,n,N,theta_g,kpq)

tic;
% n = 20;
% N = 20;

a0 = theta(1); % exact value alpha
p0 = theta(2); % exact value p
q0 = theta(3); % exact value q

% Parameter guess
ag = theta_g(1);
pg = theta_g(2);
qg = theta_g(3);

% check the condition (q>p) in the initial condition
if qg<pg 
    t=qg;
    qg=pg;
    pg=t;
end


ae = 1000;
pe = 1000;
qe = 1000;

%% EM algorithm 

% generate N data sequence of binomial mixture model

% kp = binornd(n,p0,N,1)'; % generate N samples from binomial dist with p
% kq = binornd(n,q0,N,1)'; % generate N samples from binomial dist with q
% 
% ind = [rand(N,1)'<a0];  % indicator to select either kp or kq with prob a
% kpq = ind .* kp + (1 - ind) .* kq; %mixture model at value alpha=a

e = 1e-5;                % delta between solution to stop

[row,col]=size(kpq);

maxiter=500;

next = false;

while (abs(ae-ag)>e && abs(pe-pg)>e && abs(qe-qg)>e && maxiter) 
% a0
% ae 
% ag
    if (next)
        ag = ae;
        pg = pe;
        qg = qe;
    end
    
Pl1 = [];
Pl2 = [];
for i=[1:row]
    k = kpq(i,1);
    Pp = binopdf(k,n,pg); % the first binomial distribution (p_guess)
    Pq = binopdf(k,n,qg); % the secnond binomial distribution (q_guess)

    P = ag .* Pp + (1-ag) .* Pq; % the mixture model
    
    Pl1 = [Pl1 ; (ag .* Pp)./P];         % using Baye's rule prob that observed k comes from Pp
    Pl2 = [Pl2 ; ((1-ag) * Pq)/P];       % prob that observed k comes from Pq
end

ae = sum(Pl1)/N;
pe = sum(Pl1.*kpq) / (sum(Pl1.*n));
qe = sum(Pl2.*kpq) / (sum(Pl2.*n));

next = true;

% if pe>qe
%     t=1
% else
%     t=0
% end

maxiter=maxiter-1;
end

fprintf('Converged in %d steps.\n',500-maxiter);


% check the condition (q>p) in the final step
if qe<pe 
    t=qe;
    qe=pe;
    pe=t;
    ae=1-ae;
end


theta_hat=[ae;pe;qe];
error=theta_hat-theta;

exec_time=toc;
