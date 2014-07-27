% Created by: Sherif Abdelwahab
% Updated: Nov 20, 2013
% Description: Perform intitializations

function [FIMd,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,mCRBp,mCRBq,momtime,...
    momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q]=initsim(row)

FIMd=ones(row,1)*NaN;
CRBd=ones(row,1)*NaN;
CRBa=ones(row,1)*NaN;
CRBp=ones(row,1)*NaN;
CRBq=ones(row,1)*NaN;

mCRBd=[];
mCRBa=[];
mCRBp=[];
mCRBq=[];
momtime=[];
momerror_a=[];
momerror_p=[];
momerror_q=[];
emtime=[];
emerror_a=[];
emerror_p=[];
emerror_q=[];



