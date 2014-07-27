% Created: Tamara AlShammari
% Updated: Nov 20, 2013

function [Matrix]= T(k,p,q,alpha,n)
% for fixed one alpha 
% taking just k1 >>> i.e. for 1st sequence
f_p = binopdf(k,n,p);
f_q = binopdf(k,n,q);
dalpha = (f_p - f_q)/(alpha*f_p+(1-alpha)*f_q);
dp = (alpha*f_p*((k-n*p)/(p-p^2)))/((alpha*f_p)+((1-alpha)*f_q));
dq = ((1-alpha)*f_q*((k-n*q)/(q-q^2)))/((alpha*f_p)+((1-alpha)*f_q));
Matrix=[dalpha^2 dalpha*dp dalpha*dq; dp*dalpha dp^2 dp*dq; dq*dalpha dq*dp dq^2];
end
