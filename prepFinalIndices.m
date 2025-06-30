function [tvec,f]=prepFinalIndices(tab,hosp)
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart=cumsum(monthDur);
xlabels={'Jan 20','Mar 20','May 20','Jul 20','Sep 20','Nov20','Jan 21','Mar 21','May 21','Jul 21','Sep 21','Nov21','Jan 22'};

tvec=[1;find(tab{:,5}==1);tab{end,1}+1];
difft=diff(tvec);
ld=length(difft);
inds=(1:ld)';
inds=repelem(inds,difft);

mat=table2array(tab(:,3:4));
mat=[accumarray(inds,mat(:,1),[inds(end),1],@(x)mean(x,'omitnan')),accumarray(inds,-mat(:,2),[inds(end),1],@(x)mean(x,'omitnan'))];
mat2=[repelem(mat(:,1),difft),repelem(mat(:,2),difft)];%Daily values of aggregated/repeated indices
mat=[zeros(3,2);mat(2:end,:)]';
f=mat;
tvec=[monthStart(1:3),tvec(2:end)'];

thosp=85:length(hosp)+84;
tvech=tvec(tvec<=thosp(end))';
tvech=tvech(4:end);
hval=hosp(tvech-84);
diffth=diff(tvech);
hosp=repelem(hval(1:end-1),diffth)/5e3;
thosp=85:length(hosp)+84;

mat2(1:91,:)=0;
hosp(thosp<92)=0;

fs=12; lw=2;
cmap=lines(7);
h=[];
figure
hold on
plot([0,monthStart(end)],[0,0],'k-','linewidth',1)
h(1)=plot(mat2(:,1),'-','linewidth',lw,'color',cmap(1,:));
h(2)=plot(mat2(:,2),'-','linewidth',lw,'color',cmap(2,:));
h(3)=plot(thosp,hosp,'-','linewidth',lw,'color',cmap(3,:));
%plot([0,92],[0,0],'k-','linewidth',2)
plot(92*[1,1],[-1,2],'-','linewidth',lw,'color',.5*[1,1,1])
plot(474*[1,1],[-1,2],'--','linewidth',lw,'color',.5*[1,1,1])%monthStart(20)
xticks(monthStart(1:2:end))
xticklabels(xlabels)
xtickangle(45)
set(gca,'fontsize',fs)
axis([1,monthStart(end-1)+1,0,1])
xlabel('Time')
ylabel('Index')
legend(h,'% "not improving"','Stringency','Hospital admissions/5k (data)','location','SE') %Restaurant footfall (% change)
grid on
grid minor
box on

%{
tvec=table2array(tab(:,1))';
tvec=[monthStart(1:4),tvec(2:end)];
mat=table2array(tab(:,[4,6]))';

mat(1,1:2)=mat(1,3);
mat=[zeros(2,3),mat(:,1:end-1)]/100;
f=mat;
%}