function f=plotSigma%(thetaInt)
%
x=0:.001:1;%10:1e4;%:.001:1;%**Change for x/H
L=.7;%1;%.9;
k=-12;%-.0015;%-12;
x0=.5;%6e3;%.5;
m=0;%.3;
%}
%{
x=0:100:5e4;
L=.9;
k=5e-4;
x0=15e3;
m=.4;
%}
%{
x=0:100:5e4;
L=thetaInt(1);
k=thetaInt(2);
x0=thetaInt(3);
m=thetaInt(4);
%}
fun=@(xx)(m+(1-m)*(L/(1+exp(k*(xx-x0)))));
deriv=@(xx)(-(1-m)*L*k/(1+exp(k*(xx-x0)))^2/exp(k*(xx-x0)));
y=arrayfun(fun,x);
fs=15; lw=2;
cmap=lines(7);
int=5e3;
xend=x(end);%5e4;%**Change for x/H
f=figure('Units','centimeters','Position',[0 0 12 10]);%20 16
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',15);
hold on
plot([0,xend],m+(1-m)*L*[1,1],'--','linewidth',lw,'color',.7*[1,1,1])
%plot([0,xend],m*[1,1],'--','linewidth',lw,'color',.7*[1,1,1])
plot(x0*[1,1],[0,fun(x0)],'--','linewidth',lw,'color',.7*[1,1,1])
plot(x0+int*[-1,1],fun(x0)+deriv(x0)*int*[-1,1],'linewidth',lw,'color',.7*[1,1,1])
plot(x,y,'linewidth',3,'color',cmap(1,:))
%
%
%text(.01,m+(1-m)*L-.03,'Max value $m+L(1-m)$','fontsize',fs);%.7*xend   ...L-.05 %**Change for x/H
%text(.8*xend,.33,'Min value $m$','fontsize',fs);
%text(.51,.03,'Mid point $v_0$','fontsize',fs);
%text(.25,.8,'Gradient $-kL(1-m)$','fontsize',fs);%.27*xend,.85
text(.8,m+(1-m)*L+.05,'Max. $m$','fontsize',fs);%.7*xend   ...L-.05 %**Change for x/H %.01,m+(1-m)*L-.05
%text(.65*xend,.35,'Min. value $m$','fontsize',fs);
text(.55,.05,'Mid point $v_0$','fontsize',fs);%.52,.05
text(.45,.85,'Grad. $km$','fontsize',fs);%.27*xend,.85 %.01,.65
%}
%{
text(.6*xend,m+(1-m)*L-.04,'Max value $m+L(1-m)$','fontsize',fs);%.7*xend   ...L-.05
text(.02*xend,.35,'Min value $m$','fontsize',fs);
text(.32*xend,.05,'Mid point $x_0$','fontsize',fs);
text(.4*xend,.3,'Gradient $-kL(1-m)$','fontsize',fs);%.27*xend,.85
%}
%set(gca,'fontsize',fs)
xlabel('$v(t)$')%('$H(t)$')%**Change for x/H
ylabel('$p[v(t)]$')%('$\sigma[H(t)]$')%**Change for x/H
axis ([0,x(end),0,1])%([0,5e4,0,1])%**Change for x/H
grid on
box on

