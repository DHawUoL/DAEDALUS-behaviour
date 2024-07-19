function f=prepMobilty(mobTab)
%Directly imported mobility data table
mobMat=table2array(mobTab(:,2:end));%Removes column headers
[coeff2,score2,latent2]=pca(mobMat);
xmu=mean(mobMat,1);
%plot((coeff2(:,1:2)'*(xTry-xmu))')
%xpcs=(mobMat-repmat(xmu,size(mobMat,1),1))*coeff2(:,1:3);

monthDur=[30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31];%From April 2020
months=cumsum([1,monthDur]);
lt=length(monthDur);

%1/3/20 = row 10
%31/07/2021 = row 527
%31/1/22 = row 711

x=mobMat(41:711,:);
[a,b]=size(x);
Xmonth=zeros(lt,b);
for i=1:lt
    Xmonth(i,:)=mean(x(months(i):months(i+1)-1,:),1);
end
Xmonth=[ones(b,2),Xmonth'/100];
f=Xmonth;
