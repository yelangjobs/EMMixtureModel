clear all; clc;
runs = 1000; 
n = 20;
N = 20;
 
p = 0.2;
q = 0.4;
 
alpha = [0.1:0.1:0.9]; 
%[sz1,sz2]=size(alpha);

%% Monte Carlo Simulation
for a=alpha
    Matrix=zeros(3,3);
    for r=1:runs
        kp = binornd(n,p,N,1); % generate N samples from binomial dist with p
        kq = binornd(n,q,N,1); % generate N samples from binomial dist with q
        ind = [rand(N,1)<a];  % indicator to select either kp or kq with prob a
        kpq = ind .* kp + (1 - ind) .* kq;
        Matrix=Matrix+T(kpq(1),p,q,a,n);
    end
    avgFIM=Matrix/runs
end

%% Theoritical Result
for a=alpha
    Matrix=zeros(3,3);
    for k=0:n
        f_q = binopdf(k,n,q);
        f_p = binopdf(k,n,p);
        f_pq=a * f_p + ( 1 - a ) * f_q;
        c=f_pq.*T(k,p,q,a,n);
        Matrix=Matrix+c;
    end
    tFIM=Matrix
end