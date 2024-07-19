function f=plotLasagne(x)
%%
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];
monthStart=cumsum(monthDur);
%tvec=[-68.7792,32,88.2788,monthStart(5:end)];
tvec=[0,monthStart(2),92.5467,monthStart(5:end)];
%
[a,b]=size(x);

xvec=.5:2:b+.5;%tvec(1:2:end);
xlabels={'Jan 2020','Lockdown','Jun 2020','Aug 2020','Oct 2020','Dec 2020','Feb 2021','Apr 2021','Jun 2021','Aug2021'};
%yvec=0;
fs=10; lw=2;

ff=figure;%('Units','centimeters','Position',[0 0 20 12]);
%{
set(ff,'defaulttextInterpreter','latex');
set(ff,'defaultAxesTickLabelInterpreter','latex');
set(ff,'defaultLegendInterpreter','latex');
set(ff,'DefaultAxesFontSize',fs);
%}
hold on
%{
xmin=min(min(x));
h=pcolor(1:b+1,1:a,[x,xmin*ones(a,1)]);%tc,xc,
set(h, 'EdgeColor', 'none');
%}
imagesc(x)
for j=.5:1:b+.5
    plot([j,j],[.5,a+.5],'-','linewidth',1,'color','k')
end
plot([.5,b+.5],[.5,.5],'-','linewidth',1,'color','k')
plot([.5,b+.5],[a+.5,a+.5],'-','linewidth',1,'color','k')
%}
title('GVA (fraction of pre-pandemic value)')
xlabel('Period')
ylabel('')
set(gca,'fontsize',fs)
set(gca,'xtick',xvec,'xticklabels',xlabels,'XTickLabelRotation',45)%,'ytick',yvec,'yticklabels',ylabels);
colorbar
caxis([0,1])
axis tight
box on
