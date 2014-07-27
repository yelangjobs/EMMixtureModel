%% Project Assignment 5 - Estimating Multinomial mixture models
% Created By: Sherif Abdelwahab, Tamara AlShammari
% Updated: Dec 08, 2013, Sherif

% Cleanup
clear, clc, format compact, close all


%% parameter definitions
 
runs = 200; % change to 200

n = 20;
p = 0.2; % changed this for MOM to have q>p
q = 0.4; % changed this for MOM to have q>p
N = 200; % changed this for MOM to be valid
a = 0.3;

P = [0.1:0.1:0.9]';
Q = [0.1:0.1:0.9]';
 
alpha = [0.1:0.2:0.9]';                     

ns=[3:3:30]';

Ns = [3:50:453]';

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 11)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

%% Get the CRLB in function of alpha

% theory
[row,col]=size(alpha);
%initializations
[FIMd,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,mCRBp,mCRBq,momtime,...
    momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q]=initsim(row);

for r=1:row
    ai=alpha(r,1);   
    [FIM,CRB,FIMd(r,1),CRBd(r,1),CRBa(r,1),CRBp(r,1),CRBq(r,1)]=...
        computefim(n,p,q,ai,N);
 end

% simulation
for r=1:row
  %  r
    ai=alpha(r,1);    
    Matrix=zeros(3,3);
    for r=1:runs
        kpq=mixture(n,N,p,q,ai);
        
        
        Matrix=Matrix+N*T(kpq(1),p,q,ai,n);

        
        % the mom
        [theta_hat,error,exec_time]=mom([ai;p;q],n,N,rand(3,1),kpq);
        momtime = [momtime exec_time];
        error2 = error*error';
        momerror_a = [momerror_a error2(1,1)];
        momerror_p = [momerror_p error2(2,2)];
        momerror_q = [momerror_q error2(3,3)];
        
        % the em algorithm
        [theta_hat,error,exec_time]=em([ai;p;q],n,N,rand(3,1),kpq);
        emtime = [emtime exec_time];
        error2 = error*error';
        emerror_a = [emerror_a error2(1,1)];
        emerror_p = [emerror_p error2(2,2)];
        emerror_q = [emerror_q error2(3,3)];
        
    end
    avgFIM=Matrix/runs;
    avgCRB=inv(avgFIM);
    mCRBd=[mCRBd 1/det(avgFIM)];
    mCRBa=[mCRBa avgCRB(1,1)];
    mCRBp=[mCRBp avgCRB(2,2)];
    mCRBq=[mCRBq avgCRB(3,3)];
end

Xaxis=alpha;
Xlabel='\alpha';
Simpars='$\alpha=0.1..0.9, p=0.2, q=0.4, n=20$';

plotprojfig(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'a')

plotprojfigmean(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'a')

% %% Get the CRLB in function of p
% 
% % theory
% [row,col]=size(P);
%     
% %initializations
% [FIMd,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,mCRBp,mCRBq,momtime,...
%     momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
%     emerror_p,emerror_q]=initsim(row);
% 
% for r=1:row
%     pi=P(r,1);   
%     [FIM,CRB,FIMd(r,1),CRBd(r,1),CRBa(r,1),CRBp(r,1),CRBq(r,1)]=...
%         computefim(n,pi,q,a,N);
%  end
% 
% % simulation
% for r=1:row
%     r
%     pi=P(r,1);    
%     Matrix=zeros(3,3);
%     for r=1:runs
%         kpq=mixture(n,N,pi,q,a);
%         
%         
%         Matrix=Matrix+N*T(kpq(1),pi,q,a,n);
% 
%         
%         % the mom
%         [theta_hat,error,exec_time]=mom([a;pi;q],n,N,rand(3,1),kpq);
%         momtime = [momtime exec_time];
%         error2 = error*error';
%         momerror_a = [momerror_a error2(1,1)];
%         momerror_p = [momerror_p error2(2,2)];
%         momerror_q = [momerror_q error2(3,3)];
%         
%         % the em algorithm
%         [theta_hat,error,exec_time]=em([a;pi;q],n,N,rand(3,1),kpq);
%         emtime = [emtime exec_time];
%         error2 = error*error';
%         emerror_a = [emerror_a error2(1,1)];
%         emerror_p = [emerror_p error2(2,2)];
%         emerror_q = [emerror_q error2(3,3)];
%         
%     end
%     avgFIM=Matrix/runs;
%     avgCRB=inv(avgFIM);
%     mCRBd=[mCRBd 1/det(avgFIM)];
%     mCRBa=[mCRBa avgCRB(1,1)];
%     mCRBp=[mCRBp avgCRB(2,2)];
%     mCRBq=[mCRBq avgCRB(3,3)];
% end
% 
% 
% Xaxis=P;
% Xlabel='p';
% Simpars='$p=0.1..0.9, \alpha=0.3, q=0.4, n=20$';
% 
% plotprojfig(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
%     mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
%     emerror_p,emerror_q,runs,'p')
% 

%% Get the CRLB in function of q

% theory
[row,col]=size(Q);
    

%initializations
[FIMd,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,mCRBp,mCRBq,momtime,...
    momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q]=initsim(row);

for r=1:row
    qi=Q(r,1);   
    [FIM,CRB,FIMd(r,1),CRBd(r,1),CRBa(r,1),CRBp(r,1),CRBq(r,1)]=...
        computefim(n,p,qi,a,N);
end

%simulation

for r=1:row
    r
    qi=Q(r,1);    
    Matrix=zeros(3,3);
    for r=1:runs
        kpq=mixture(n,N,p,qi,a);
        Matrix=Matrix+N*T(kpq(1),p,qi,a,n);
        
        
        % the mom
        [theta_hat,error,exec_time]=mom([a;p;qi],n,N,rand(3,1),kpq);
        momtime = [momtime exec_time];
        error2 = error*error';
        momerror_a = [momerror_a error2(1,1)];
        momerror_p = [momerror_p error2(2,2)];
        momerror_q = [momerror_q error2(3,3)];
        
        % the em algorithm
        [theta_hat,error,exec_time]=em([a;p;qi],n,N,rand(3,1),kpq);
        emtime = [emtime exec_time];
        error2 = error*error';
        emerror_a = [emerror_a error2(1,1)];
        emerror_p = [emerror_p error2(2,2)];
        emerror_q = [emerror_q error2(3,3)];
        
        
    end
    avgFIM=Matrix/runs;
    avgCRB=inv(avgFIM);
    mCRBd=[mCRBd 1/det(avgFIM)];
    mCRBa=[mCRBa avgCRB(1,1)];
    mCRBp=[mCRBp avgCRB(2,2)];
    mCRBq=[mCRBq avgCRB(3,3)];
end

Xaxis=Q;
Xlabel='q';
Simpars='$q=0.1..0.9, \alpha=0.3, p=0.2, n=20$';

plotprojfig(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'q')


plotprojfigmean(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'q')



%% Get the CRLB in function of n

%theory
[row,col]=size(ns);
    

%initializations
[FIMd,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,mCRBp,mCRBq,momtime,...
    momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q]=initsim(row);
    
for r=1:row
    ni=ns(r,1);   
    [FIM,CRB,FIMd(r,1),CRBd(r,1),CRBa(r,1),CRBp(r,1),CRBq(r,1)]=...
        computefim(ni,p,q,a,N);
end

%simulation
for r=1:row
    r
    ni=ns(r,1);    
    Matrix=zeros(3,3);
    for r=1:runs
        kpq=mixture(ni,N,p,q,a);
        Matrix=Matrix+N*T(kpq(1),p,q,a,ni);
        
        
        % the mom
        [theta_hat,error,exec_time]=mom([a;p;q],ni,N,rand(3,1),kpq);
        momtime = [momtime exec_time];
        error2 = error*error';
        momerror_a = [momerror_a error2(1,1)];
        momerror_p = [momerror_p error2(2,2)];
        momerror_q = [momerror_q error2(3,3)];
        
        % the em algorithm
        [theta_hat,error,exec_time]=em([a;p;q],ni,N,rand(3,1),kpq);
        emtime = [emtime exec_time];
        error2 = error*error';
        emerror_a = [emerror_a error2(1,1)];
        emerror_p = [emerror_p error2(2,2)];
        emerror_q = [emerror_q error2(3,3)];
        
            
        
    end
    avgFIM=Matrix/runs;
    avgCRB=inv(avgFIM);
    mCRBd=[mCRBd 1/det(avgFIM)];
    mCRBa=[mCRBa avgCRB(1,1)];
    mCRBp=[mCRBp avgCRB(2,2)];
    mCRBq=[mCRBq avgCRB(3,3)];
end



Xaxis=ns;
Xlabel='n';
Simpars='$n=3..30, \alpha=0.3, p=0.2, q=0.4$';

plotprojfig(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'n')



plotprojfigmean(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'n')


%% Get the CRLB in function of N

%theory
[row,col]=size(Ns);
    

%initializations
[FIMd,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,mCRBp,mCRBq,momtime,...
    momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q]=initsim(row);
    
for r=1:row
    Ni=Ns(r,1);   
    [FIM,CRB,FIMd(r,1),CRBd(r,1),CRBa(r,1),CRBp(r,1),CRBq(r,1)]=...
        computefim(n,p,q,a,Ni);
end

a_hat=[];
%simulation
for r=1:row
    r
    Ni=Ns(r,1);    
    Matrix=zeros(3,3);
    myas=0;
    for r=1:runs
        kpq=mixture(n,Ni,p,q,a);
        Matrix=Matrix+Ni*T(kpq(1),p,q,a,n);
        
        
        % the mom
        [theta_hat,error,exec_time]=mom([a;p;q],n,Ni,rand(3,1),kpq);
        if theta_hat(1)==0 
            myas=myas+1;
        end
        momtime = [momtime exec_time];
        error2 = error*error';
        momerror_a = [momerror_a error2(1,1)];
        momerror_p = [momerror_p error2(2,2)];
        momerror_q = [momerror_q error2(3,3)];
        
        % the em algorithm
        [theta_hat,error,exec_time]=em([a;p;q],n,Ni,rand(3,1),kpq);
        emtime = [emtime exec_time];
        error2 = error*error';
        emerror_a = [emerror_a error2(1,1)];
        emerror_p = [emerror_p error2(2,2)];
        emerror_q = [emerror_q error2(3,3)];
        
            
        
    end
    avgFIM=Matrix/runs;
    avgCRB=inv(avgFIM);
    mCRBd=[mCRBd 1/det(avgFIM)];
    mCRBa=[mCRBa avgCRB(1,1)];
    mCRBp=[mCRBp avgCRB(2,2)];
    mCRBq=[mCRBq avgCRB(3,3)];
    a_hat=[a_hat myas/runs];
end

Xaxis=Ns;
Xlabel='N';
Simpars='$N=10..100, \alpha=0.3, p=0.2, n=20$';


h=figure()
plot(Ns,a_hat); grid on;
xlabel(Xlabel)
ylabel('Indeterminant estimator Prob','Interpreter','LaTex')
%title(strcat('Performance of MOM estimator of $\hat{q}$ vs ',Simpars),'interpreter','latex')
saveas(h,strcat('momproblem_','Num'),'fig')
saveas(h,strcat('momproblem_','Num'),'eps')



plotprojfig(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'Num')


plotprojfigmean(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,'Num')


