tic;
n = 20;
N = 20;

a0 = 0.8; % exact value alpha
p0 = 0.2; % exact value p
q0 = 0.4; % exact value q




%% MOM algorithm 

% generate N data sequence of binomial mixture model

kp = binornd(n,p0,N,1)'; % generate N samples from binomial dist with p
kq = binornd(n,q0,N,1)'; % generate N samples from binomial dist with q

ind = [rand(N,1)'<a0];  % indicator to select either kp or kq with prob a
kpq = ind .* kp + (1 - ind) .* kq; %mixture model at value alpha=a



[row,col]=size(kpq);

M1 = [];
M2 = [];
M3 = [];

for i=[1:col]
    k = kpq(1,i);
     
    M1 = [M1 k ];
    M2 = [M2 k * (k-1) ];
    M3 = [M3 k * (k-1) * (k-2) ];
       
    
end

M1 = sum(M1) / N / n;
M2 = sum(M2) / N / n / (n-1);
M3 = sum(M3) / N / n / (n-1) / (n-3);

Z = (M3 - M2 * M1) / (M2 - M1^2)

root = sqrt(Z^2-4*M1*Z+4*M2)

pe = 0.5 * (Z - root)

qe = 0.5 * (Z + root)

ae = (M1 - qe)/(pe - qe)


toc