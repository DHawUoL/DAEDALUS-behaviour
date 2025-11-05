function f=bePlotIndices(Xdaily,tvec,hosp)
hospStart=85;%Start day of admisisons data
Xdaily(1,1:91)=0;
figure
fs=12; lw=2;
hold on
plot(Xdaily','linewidth',2);
plot(85:tvec(end),hosp(1:tvec(end)-hospStart+1)/5e3,'linewidth',lw);
legend('Trust','Stringency','Hospital admissions/4k','location','SW')
set(gca,'fontsize',fs)
box on
grid on
grid minor