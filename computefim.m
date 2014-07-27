% Created by: Sherif Abdelwahab
% Updated: Nov 20, 2013
% Compute the theoritical FIM and CRB

function [FIM, CRB, FIMd,CRBd,CRBa,CRBp,CRBq] = computefim(n,p,q,a,N)

FIMd=NaN;
CRBd=NaN;
CRBa=NaN;
CRBp=NaN;
CRBq=NaN;

if p==q
    disp('FIM is not PSD');
    FIM=zeros(3,3);
    CRB=eye(3);
    return
end

c=[];
for k=[0:n]
    c=[c nchoosek(n,k)];

end
k=[0:n];

Pp =  c .* power(p,k) .* power((1-p),(n-k));

Pq =  c .* power(q,k) .* power((1-q),(n-k));

P = a .* Pp + (1-a) .* Pq;

Da = power(P,-1) .* (Pp - Pq);

Dp = power(P,-1) .* a .* Pp .* ((k-n*p)/p/(1-p));

Dq = power(P,-1) .* (1-a) .* Pq .* ((k-n*q)/q/(1-q));

Eaa = sum(Da .* Da .* P);
Eap = sum(Da .* Dp .* P);
Eaq = sum(Da .* Dq .* P);
Epp = sum(Dp .* Dp .* P);
Epq = sum(Dp .* Dq .* P);
Eqq = sum(Dq .* Dq .* P);

FIM = N*[ Eaa Eap Eaq; Eap Epp Epq; Eaq Epq Eqq];
CRB = inv(FIM);


FIMd=det(FIM);
    if p~=q   
        CRBd=det(CRB);
        CRBa=CRB(1,1);
        CRBp=CRB(2,2);
        CRBq=CRB(3,3);
    end
    
    
    