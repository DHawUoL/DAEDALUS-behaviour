function f=expandStringency(x,data)
%x: day number, value
[a,b]=size(x);
xend=x(a,1);
xout=cumsum(ones(xend,1));

%x1:
addCol=zeros(xend,1);
for i=1:a-1
    addCol(x(i,1):x(i+1,1)-1)=x(i,2);%b 2
end
xout=[xout,addCol];
addCol=zeros(xend,1);
%x2:
addCol=zeros(xend,1);
for i=1:a-1
    addCol(x(i,1):x(i+1,1)-1)=x(i,3);
end
xout=[xout,addCol];
addCol=zeros(xend,1);
%x3:
addCol=zeros(xend,1);
for i=1:a-1
    addCol(x(i,1):x(i+1,1)-1)=x(i,4);
end
xout=[xout,addCol];

f=xout;
fs=12; lw=2; 
fromHere=85;
dataTop=1e3*ceil(max(data)/1e3);
plotend=min(xend,length(data)+fromHere-1);

figure
hold on
bar(fromHere:plotend,data(1:plotend-fromHere+1)/dataTop)
plot(xout(fromHere:plotend,1)+14,xout(fromHere:plotend,2:end)/100,'linewidth',lw)
set(gca,'fontsize',fs)
xlabel('Day')
ylabel('Variable')
grid on
box on
