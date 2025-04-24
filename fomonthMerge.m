function f=fomonthMerge(tab)%mobility table input
monthDur=[1,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31];
monthStart=cumsum(monthDur);

numDays=[repelem((1:length(monthStart)-1)',monthDur(2:end)');length(monthStart)];
numDays=numDays(tab{1,3}:tab{end,3});
%x=[numDays,tab{:,4}+100];
x=accumarray(numDays,tab{:,4}+100,[],@mean);
f=x;