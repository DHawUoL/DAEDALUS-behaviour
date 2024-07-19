function f=plotLasagneIndices(dat)
xlab={'Jan 20','Mar 20','May 20','Jul 20','Sep 20','Nov 20','Jan 21','Mar 21','May 21','July 21'};
ylab={'Stringency','Government Response','Containment and health',...
    'Economic support','School closures','Workplace closures'...
    'Cancelling public events','Limits on gatherings','Closing of public transport',...
    'Stay-at-home requirements','Restrictions on internal movement','Restrictions on international travel'};
[a,b]=size(dat);
imagesc(dat');
xticks(1:2:a)
xticklabels(xlab)
xtickangle(45)
yticks(1:b)
yticklabels(ylab)