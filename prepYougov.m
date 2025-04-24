function f=prepYougov(tab)
%Dates/indices - weekly yougov data, confidence

%Dates to numeric (
tab.day=convertTo(tab{:,1},'datenum')-convertTo(datetime('31-Dec-2019'),'datenum');

%Manually fill in missing line:
tab(10,2:end-1)=tab(9,2:end-1);

fs=10; lw=2; ms=2;
cmap=lines(7);
figure
colormap(cmap);
plot(tab.day,tab{:,2:end-1},'-o','linewidth',lw)%'markersize',ms,'markerfacecolor',cmap(1,:))
%rod(tab{:,2:end-1},2)

%xlim([xdata(1),xdata(end)]);
%ylim([0,ymax])%.25*max(ydata)]);
%xticks(xvec)
%xticklabels(xnames)
%xtickangle(45)
box on;
grid on;
xlabel('Week');
ylabel('Index');
%title('Model Fit');