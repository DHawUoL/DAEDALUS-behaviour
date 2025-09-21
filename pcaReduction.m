function f=pcaReduction(x)
%xfull=x(:,4:60)'; %Non-zero

tvec=[1    32    61    93   107   169   176   223   227   230   250   258   265   271   279   294   310   322   330   338   349   370   384    407   418   433   445   463   474   491   504   517   540   561   567   575   594   605   617   631   642   652   661   679   693   702    716   742   784];
lt=length(tvec);
tvec=tvec(1:(lt-15));
tdiff=diff(tvec);

xfull=x(:,4:length(tdiff))';%4:60
%xfull=xfull(1:length(tdiff),:);
xdiff=[diff(xfull(:,1)),diff(xfull(:,2))];
exclude=sum(abs(xdiff),2);
exclude=find(exclude==0);
xin=xfull;
xin(exclude+1,:)=NaN;

[coeff,score,latent,tsquared,explained,mu] = pca(xin);

lx=length(x);
xfull1=(x'-repmat(mu,lx,1))*coeff(:,1);
xfull2=(x'-repmat(mu,lx,1))*coeff(:,2);
xfull=[xfull1,xfull2]';
f=xfull;

%
%% Repeating indices:
xplot=[];
xorig=[];
for i=1:length(tdiff)
    xplot=[xplot,repmat(xfull(:,i),1,tdiff(i))];
    xorig=[xorig,repmat(x(:,i),1,tdiff(i))];
end
figure
lw=2;
hold on
plot(xplot','-','LineWidth',lw)
plot(xorig','--','LineWidth',lw)
legend('PC1','PC2','Trust','Stringency','location','northeastoutside')
%axis([1,500,0,1])
box on
grid on
%}
end