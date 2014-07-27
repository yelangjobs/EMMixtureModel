% Created By: Sherif Abdelwahab
% Updated: Dec 09, 2013
% DEscription: Plot figures with mean

function plotprojfigmean(Xaxis,Xlabel,Simpars,CRBd,CRBa,CRBp,CRBq,mCRBd,mCRBa,...
    mCRBp,mCRBq,momtime,momerror_a,momerror_p,momerror_q,emtime,emerror_a,...
    emerror_p,emerror_q,runs,suffex)

maxerror=0.5;

pl = 100;

h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1)
% plot (Xaxis, mCRBd,Xaxis,CRBd,'LineWidth',2)
% 
% grid on
% xlabel(Xlabel)
% ylabel('det(\bf{CRLB})')
% title(strcat('Determinant of Cramer-Rao lower bound vs ',Simpars),'interpreter','latex')
% legend('det(CRLB) - Monte-Carlo','det(CRLB) - Exact')    
% ylim([min([mCRBd CRBd']) max([mCRBd CRBd'])])
% 
% subplot(2,1,2)
plot(Xaxis, CRBa,'b', Xaxis, CRBp,'r',Xaxis, CRBq,'g',Xaxis, mCRBa,'b--', Xaxis, mCRBp,'r--',Xaxis, mCRBq,'g--','LineWidth',2)

grid off
xlabel(Xlabel)
ylabel('$\bf{CRLB}(\hat{\theta})$','Interpreter','LaTex')
title(strcat('Cramer-Rao lower bound of each parameter vs ',Simpars),'Interpreter','LaTex')
legend('Exact CRLB(\alpha)','Exact CRLB(p)','Exact CRLB(q)','Monte-carlo CRLB(\alpha)','Monte-carlo CRLB(p)','Monte-carlo CRLB(q)')
set(gca,'yscale','log')
saveas(h,strcat('meanCRB_',suffex),'fig')
saveas(h,strcat('meanCRB_',suffex),'eps')


h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%subplot(3,2,1)
errordata=reshape(emerror_a,runs,length(Xaxis));
boxplot(errordata,Xaxis,'extrememode','compress','datalim',[0 max(prctile(errordata,pl))]);
hold on; plot(CRBa,'g','LineWidth',2);
plot(mean(errordata),'r','LineWidth',2); grid off
xlabel(Xlabel)
ylabel('$MSE(\hat{\alpha})$','Interpreter','LaTex')
title(strcat('Performance of EM estimator of $\hat{\alpha}$ vs ',Simpars),'Interpreter','LaTex')
%legend('CRLB(\alpha)','EM MSE(\alpha)')
set(gca,'yscale','log')
ylim([0 max([prctile(errordata,pl)])])


saveas(h,strcat('meanEMvsMOM_1_',suffex),'fig')
saveas(h,strcat('meanEMvsMOM_1_',suffex),'eps')

h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%subplot(3,2,2)
errordata=reshape(momerror_a,runs,length(Xaxis));
boxplot(errordata,Xaxis,'extrememode','compress','datalim',[0 max(prctile(errordata,pl))]);
hold on; plot(CRBa,'g','LineWidth',2);
plot(mean(errordata),'r','LineWidth',2); grid off
xlabel(Xlabel)
ylabel('$MSE(\hat{\alpha})$','Interpreter','LaTex')
title(strcat('Performance of MOM estimator of $\hat{\alpha}$ vs ',Simpars),'Interpreter','LaTex')
legend('$CRLB(\alpha)$','$MOM MSE(\alpha)$','Interpreter','latex')
ylim([0 max([prctile(errordata,pl)])])
set(gca,'yscale','log')
h = legend;
set(h, 'interpreter', 'latex')


saveas(h,strcat('meanEMvsMOM_2_',suffex),'fig')
saveas(h,strcat('meanEMvsMOM_2_',suffex),'eps')

h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%subplot(3,2,3)
errordata=reshape(emerror_p,runs,length(Xaxis));
boxplot(errordata,Xaxis,'extrememode','compress','datalim',[0 max(prctile(errordata,pl))]);
hold on; plot(CRBp,'g','LineWidth',2);
plot(mean(errordata),'r','LineWidth',2); grid off
xlabel(Xlabel)
ylabel('$MSE(\hat{p})$','Interpreter','LaTex')
title(strcat('Performance of EM estimator of $\hat{p}$ vs ',Simpars),'Interpreter','LaTex')
%legend('CRLB(p)','EM MSE(\hat{p})')
ylim([0 max([prctile(errordata,pl)])])
set(gca,'yscale','log')


saveas(h,strcat('meanEMvsMOM_3_',suffex),'fig')
saveas(h,strcat('meanEMvsMOM_3_',suffex),'eps')

h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%subplot(3,2,4)
errordata=reshape(momerror_p,runs,length(Xaxis));
boxplot(errordata,Xaxis,'extrememode','compress','datalim',[0 max(prctile(errordata,pl))]);
hold on; plot(CRBp,'g','LineWidth',2);
plot(mean(errordata),'r','LineWidth',2); grid off
xlabel(Xlabel)
ylabel('$MSE(\hat{p})$','Interpreter','LaTex')
title(strcat('Performance of MOM estimator of $\hat{p}$ vs ',Simpars),'Interpreter','LaTex')
legend('$CRLB(p)$','MOM $MSE(\hat{p})$','Interpreter','LaTex')
ylim([0 max([prctile(errordata,pl)])])
h = legend;
set(h, 'interpreter', 'latex')
set(gca,'yscale','log')


saveas(h,strcat('meanEMvsMOM_4_',suffex),'fig')
saveas(h,strcat('meanEMvsMOM_4_',suffex),'eps')

h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%subplot(3,2,5)
errordata=reshape(emerror_q,runs,length(Xaxis));
boxplot(errordata,Xaxis,'extrememode','compress','datalim',[0 max(prctile(errordata,pl))]);
hold on; plot(CRBq,'g','LineWidth',2);
plot(mean(errordata),'r','LineWidth',2); grid off
xlabel(Xlabel)
ylabel('$MSE(\hat{q})$','Interpreter','LaTex')
title(strcat('Performance of EM estimator of $\hat{q}$ vs ',Simpars),'Interpreter','LaTex')
%legend('CRLB(q)','EM MSE(\hat{q})')
ylim([0 max([prctile(errordata,pl)])])
set(gca,'yscale','log')


saveas(h,strcat('meanEMvsMOM_5_',suffex),'fig')
saveas(h,strcat('meanEMvsMOM_5_',suffex),'eps')


h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%subplot(3,2,6)
errordata=reshape(momerror_q,runs,length(Xaxis));
boxplot(errordata,Xaxis,'extrememode','compress','datalim',[0 max(prctile(errordata,pl))]);
hold on; plot(CRBq,'g','LineWidth',2);
plot(mean(errordata),'r','LineWidth',2); grid off
xlabel(Xlabel)
ylabel('$MSE(\hat{q})$','Interpreter','LaTex')
title(strcat('Performance of MOM estimator of $\hat{q}$ vs ',Simpars),'Interpreter','LaTex')
legend('$CRLB(q)$','MOM $MSE(\hat{q})$','Interpreter','LaTex')
ylim([0 max([prctile(errordata,pl)])])
h = legend;
set(h, 'interpreter', 'latex')
set(gca,'yscale','log')

saveas(h,strcat('meanEMvsMOM_6_',suffex),'fig')
saveas(h,strcat('meanEMvsMOM_6_',suffex),'eps')


h=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
emtimedata=reshape(emtime,runs,length(Xaxis));
momtimedata=reshape(momtime,runs,length(Xaxis));
f1 = @(x, y) bar(x, y, 0.6);
f2 = @(x, y) bar(x, y, 0.2); % even narrower
[ax,h1,h2]=plotyy(Xaxis,mean(emtimedata)',Xaxis,mean(momtimedata)',f1,f2)
set(h1,'FaceColor','r')
set(h2,'FaceColor','b')
grid off
set(get(ax(1),'Ylabel'),'String','EM runtime (s)') 
set(get(ax(2),'Ylabel'),'String','MOM runtime (s)') 
title(strcat('Comparison of MOM and EM runtime at ',Simpars),'Interpreter','LaTex')
legend('EM mean runtime','MOM mean runtime')
h = legend;
set(h, 'interpreter', 'latex')
xlabel(Xlabel)
saveas(h,strcat('meanruntime_',suffex),'fig')
saveas(h,strcat('meanruntime_',suffex),'eps')
