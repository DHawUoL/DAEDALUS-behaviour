function f=prepData(x)
age=1;
%1 row per timepoint, 1 column per weekly variable
[a,b]=size(x);

if age==1
    x=reshape(x,a/4,b*4);%Each question separated by age
end
x=x';
[a,b]=size(x);

b2=ceil(b/4);%Approx 4 weeks=1 month
xx=zeros(a,b2);
for i=1:b2-1
    xx(:,i)=mean(x(:,(4*(i-1)+1):4*i),2,'omitnan');
end
xx(:,end)=mean(x(:,4*(b2-1)+1:end),2,'omitnan');

f=xx;