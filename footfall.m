function f=footfall(tab)
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthDur=monthDur(4:26)';
inds=repelem((3:25)',monthDur);
linds=25;

tab=tab(:,[2,3,7]);
tab=unstack(tab,'roll','activity');
tab=tab(4:705,:);%From March 2020 to jan 2022

mat=table2array(tab(:,2:end));
matMonth=[accumarray(inds,mat(:,1),[linds,1],@(x)mean(x,'omitnan')),accumarray(inds,mat(:,2),[linds,1],@(x)mean(x,'omitnan')),accumarray(inds,mat(:,3),[linds,1],@(x)mean(x,'omitnan')),accumarray(inds,mat(:,4),[linds,1],@(x)mean(x,'omitnan'))];

%groupsummary(mat(:,1),inds,'nanmean')
%grpstats(tab(:,2:end),inds,'mean');
matOut=1-matMonth(1:25,:)/100;
matOut=[zeros(3,size(matOut,2));matOut(4:end,:)];
f=matOut(1:25,:);
%xinabc=[dataMat(:,[47])/100,dataMat(:,[57])]';